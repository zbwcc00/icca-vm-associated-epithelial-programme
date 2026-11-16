from __future__ import annotations

import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:\第六篇大论文")
REVISION = ROOT / "revision_round1_reviewer3_20261116"
OUT = REVISION / "analysis_01_corrected_scoring_and_anchor_ablation"
OUT.mkdir(parents=True, exist_ok=True)

ANCHORS = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
ENDOTHELIAL_CONTROLS = ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "EGFL7", "PLVAP", "CA4"]
SEED = 20261116
MODULE_SEED = 20260830
UCELL_MAX_RANK = 1500

COHORTS = {
    "GSE138709": {
        "h5ad": ROOT / "results/11_gse138709_strict_qc_map/gse138709_strict_qc_map.h5ad",
        "old_scores": ROOT / "results/15_gse138709_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "patient": "patient_id",
    },
    "GSE125449": {
        "h5ad": ROOT / "results/43_gse125449_icca_malignant_epithelial_definition/gse125449_icca_operational_malignant_epithelium.h5ad",
        "old_scores": ROOT / "results/46_gse125449_icca_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "patient": "patient_id",
    },
}


def restore_gse138709_anchor_columns(obj: ad.AnnData) -> tuple[ad.AnnData, pd.DataFrame]:
    missing = [gene for gene in ANCHORS if gene not in obj.var_names]
    if not missing:
        return obj, pd.DataFrame()
    raw_root = ROOT / "results/01_gse138709_import_v2"
    restored = {gene: np.zeros(obj.n_obs, dtype=np.float32) for gene in missing}
    obs_position = pd.Series(np.arange(obj.n_obs), index=obj.obs_names.astype(str))
    audit: list[dict[str, object]] = []
    for sample_dir in sorted(path for path in raw_root.iterdir() if path.is_dir()):
        qc_path = sample_dir / "cell_qc_unfiltered.tsv"
        gene_path = sample_dir / "genes.tsv"
        matrix_path = sample_dir / "counts_cells_by_genes.npz"
        if not (qc_path.exists() and gene_path.exists() and matrix_path.exists()):
            continue
        qc = pd.read_csv(qc_path, sep="\t")
        genes = pd.read_csv(gene_path, sep="\t").iloc[:, 0].astype(str)
        matrix = sparse.load_npz(matrix_path).tocsr()
        common = qc["cell_id"].astype(str).isin(obs_position.index)
        if not common.any():
            continue
        raw_rows = np.flatnonzero(common.to_numpy())
        target_rows = obs_position.loc[qc.loc[common, "cell_id"].astype(str)].to_numpy(int)
        gene_index = pd.Series(np.arange(len(genes)), index=genes)
        for gene in missing:
            if gene not in gene_index.index:
                continue
            values = np.asarray(matrix[raw_rows, int(gene_index[gene])].todense()).ravel().astype(np.float32)
            restored[gene][target_rows] = values
            audit.append(
                {
                    "sample_id": sample_dir.name,
                    "gene": gene,
                    "n_h5ad_cells_matched": int(len(target_rows)),
                    "n_nonzero_cells": int((values > 0).sum()),
                    "total_UMI": int(values.sum()),
                    "source": str(matrix_path.relative_to(ROOT)),
                }
            )
    counts = obj.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    counts = counts.tocsr(copy=True)
    counts.eliminate_zeros()
    extra_counts = sparse.csr_matrix(np.column_stack([restored[gene] for gene in missing]))
    all_counts = sparse.hstack([counts, extra_counts], format="csr")
    totals = np.asarray(all_counts.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    normalized = sparse.diags(1e4 / totals).dot(all_counts)
    normalized.data = np.log1p(normalized.data)
    var = pd.DataFrame(index=list(obj.var_names.astype(str)) + missing)
    augmented = ad.AnnData(X=normalized.tocsr(), obs=obj.obs.copy(), var=var)
    augmented.layers["counts"] = all_counts
    return augmented, pd.DataFrame(audit)


def select_operational_cells(cohort: str, obj: ad.AnnData, old_scores: pd.DataFrame) -> ad.AnnData:
    cell_ids = old_scores["cell_id"].astype(str)
    absent = cell_ids[~cell_ids.isin(obj.obs_names.astype(str))]
    if len(absent):
        raise ValueError(f"{cohort}: {len(absent)} old score cells are absent from the h5ad")
    return obj[cell_ids.to_numpy()].copy()


def signature_ranks(
    counts: sparse.csr_matrix,
    target_indices: list[int],
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return UCell average-tie ranks and reproducible AUCell random-tie ranks.

    Only ranks for target genes are materialised. AUCell tie positions are sampled
    without replacement inside each expression-value tie group, which is
    distributionally equivalent to randomly shuffling the complete tie group.
    """
    matrix = counts.tocsr(copy=True)
    matrix.eliminate_zeros()
    n_cells, n_genes = matrix.shape
    target_indices_array = np.asarray(target_indices, dtype=int)
    ucell_ranks = np.empty((n_cells, len(target_indices)), dtype=np.float64)
    aucell_ranks = np.empty((n_cells, len(target_indices)), dtype=np.int32)
    rng = np.random.default_rng(seed)
    for row in range(n_cells):
        start, stop = matrix.indptr[row], matrix.indptr[row + 1]
        indices = matrix.indices[start:stop]
        values = matrix.data[start:stop]
        value_lookup = dict(zip(indices.tolist(), values.tolist()))
        target_values = np.array([value_lookup.get(index, 0.0) for index in target_indices_array])
        for value in np.unique(target_values):
            target_positions = np.flatnonzero(target_values == value)
            if value == 0:
                greater = len(values)
                equal = n_genes - len(values)
            else:
                greater = int(np.count_nonzero(values > value))
                equal = int(np.count_nonzero(values == value))
            average_rank = greater + (equal + 1) / 2.0
            ucell_ranks[row, target_positions] = average_rank
            sampled_positions = rng.choice(equal, size=len(target_positions), replace=False) + 1
            aucell_ranks[row, target_positions] = greater + sampled_positions
    return ucell_ranks, aucell_ranks


def ucell_score(ranks: np.ndarray, signature_positions: list[int], max_rank: int) -> np.ndarray:
    selected = ranks[:, signature_positions].copy()
    selected[selected >= max_rank] = max_rank
    n_signature = selected.shape[1]
    rank_sum_min = n_signature * (n_signature + 1) / 2
    denominator = n_signature * max_rank - rank_sum_min
    scores = 1 - (selected.sum(axis=1) - rank_sum_min) / denominator
    scores[np.all(selected >= max_rank, axis=1)] = 0.0
    return scores


def aucell_score(ranks: np.ndarray, signature_positions: list[int], auc_threshold: int) -> np.ndarray:
    selected = ranks[:, signature_positions]
    n_signature = selected.shape[1]
    ideal = np.arange(1, n_signature + 1)
    ideal = ideal[ideal < auc_threshold]
    max_auc = np.sum(np.diff(np.r_[ideal, auc_threshold]) * np.arange(1, len(ideal) + 1))
    scores = np.zeros(selected.shape[0], dtype=float)
    for row in range(selected.shape[0]):
        x = np.sort(selected[row][selected[row] < auc_threshold])
        if len(x):
            y = np.arange(1, len(x) + 1)
            scores[row] = np.sum(np.diff(np.r_[x, auc_threshold]) * y) / max_auc
    return scores


def module_scores(obj: ad.AnnData, signatures: dict[str, list[str]]) -> pd.DataFrame:
    output = pd.DataFrame(index=obj.obs_names.astype(str))
    work = obj.copy()
    for name, genes in signatures.items():
        field = f"module_{name}"
        sc.tl.score_genes(
            work,
            gene_list=genes,
            score_name=field,
            random_state=MODULE_SEED,
            ctrl_size=50,
            n_bins=25,
        )
        output[field] = work.obs[field].astype(float).to_numpy()
    return output


def call_high(
    score_frame: pd.DataFrame,
    signature_name: str,
    vm_quantile: float,
    endothelial_quantile: float = 0.90,
) -> pd.Series:
    vm_fields = [f"module_{signature_name}", f"ucell_{signature_name}", f"aucell_{signature_name}"]
    control_fields = ["module_endothelial", "ucell_endothelial", "aucell_endothelial"]
    vm_high = pd.concat(
        [score_frame[field].gt(score_frame[field].quantile(vm_quantile)) for field in vm_fields], axis=1
    ).all(axis=1)
    endothelial_high = pd.concat(
        [score_frame[field].gt(score_frame[field].quantile(endothelial_quantile)) for field in control_fields], axis=1
    ).any(axis=1)
    return vm_high & ~endothelial_high


def kappa_binary(reference: np.ndarray, alternative: np.ndarray) -> float:
    reference = reference.astype(bool)
    alternative = alternative.astype(bool)
    observed = np.mean(reference == alternative)
    p_ref = np.mean(reference)
    p_alt = np.mean(alternative)
    expected = p_ref * p_alt + (1 - p_ref) * (1 - p_alt)
    return (observed - expected) / (1 - expected) if expected < 1 else np.nan


def classification_summary(
    cohort: str,
    labels: dict[str, pd.Series],
    patients: pd.Series,
) -> pd.DataFrame:
    full = labels["full6"].to_numpy(bool)
    rows: list[dict[str, object]] = []
    for name, label in labels.items():
        values = label.to_numpy(bool)
        union = np.logical_or(full, values).sum()
        intersection = np.logical_and(full, values).sum()
        rows.append(
            {
                "cohort": cohort,
                "classification": name,
                "omitted_anchor": name.replace("drop_", "") if name.startswith("drop_") else "none",
                "n_cells": len(values),
                "n_high": int(values.sum()),
                "high_fraction": float(values.mean()),
                "n_patients": int(patients.nunique()),
                "n_patients_with_high": int(patients[values].nunique()),
                "jaccard_vs_corrected_full6": intersection / union if union else 1.0,
                "cohen_kappa_vs_corrected_full6": kappa_binary(full, values),
                "n_gained_vs_corrected_full6": int(np.logical_and(~full, values).sum()),
                "n_lost_vs_corrected_full6": int(np.logical_and(full, ~values).sum()),
            }
        )
    return pd.DataFrame(rows)


def held_out_pseudobulk(
    cohort: str,
    obj: ad.AnnData,
    labels: dict[str, pd.Series],
    patients: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = obj.layers["counts"].tocsr()
    gene_index = pd.Series(np.arange(obj.n_vars), index=obj.var_names.astype(str))
    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for gene in ANCHORS:
        label = labels[f"drop_{gene}"].to_numpy(bool)
        column = np.asarray(counts[:, int(gene_index[gene])].todense()).ravel()
        patient_differences = []
        for patient in sorted(patients.astype(str).unique()):
            patient_mask = patients.astype(str).eq(patient).to_numpy()
            group_values = {}
            for group, group_mask in [("high", label), ("low", ~label)]:
                mask = patient_mask & group_mask
                total_umi = float(counts[mask].sum()) if mask.any() else 0.0
                gene_umi = float(column[mask].sum()) if mask.any() else 0.0
                log_cpm = math.log1p(gene_umi / max(total_umi, 1.0) * 1e6) if mask.any() else np.nan
                group_values[group] = log_cpm
                detail_rows.append(
                    {
                        "cohort": cohort,
                        "omitted_anchor": gene,
                        "patient_id": patient,
                        "classification_group": group,
                        "n_cells": int(mask.sum()),
                        "held_out_gene_UMI": gene_umi,
                        "total_UMI": total_umi,
                        "held_out_gene_log1p_CPM": log_cpm,
                    }
                )
            if np.isfinite(group_values["high"]) and np.isfinite(group_values["low"]):
                patient_differences.append(group_values["high"] - group_values["low"])
        differences = np.asarray(patient_differences, dtype=float)
        if len(differences) >= 2 and np.any(differences != 0):
            test = stats.wilcoxon(differences, alternative="two-sided", method="auto")
            p_value = float(test.pvalue)
            statistic = float(test.statistic)
        else:
            p_value = np.nan
            statistic = np.nan
        summary_rows.append(
            {
                "cohort": cohort,
                "omitted_anchor": gene,
                "n_paired_patients": len(differences),
                "mean_high_minus_low_log1p_CPM": float(np.mean(differences)) if len(differences) else np.nan,
                "median_high_minus_low_log1p_CPM": float(np.median(differences)) if len(differences) else np.nan,
                "wilcoxon_statistic": statistic,
                "nominal_p": p_value,
            }
        )
    summary = pd.DataFrame(summary_rows)
    valid = summary["nominal_p"].notna()
    summary["BH_within_cohort_6_tests"] = np.nan
    if valid.any():
        summary.loc[valid, "BH_within_cohort_6_tests"] = multipletests(
            summary.loc[valid, "nominal_p"], method="fdr_bh"
        )[1]
    return pd.DataFrame(detail_rows), summary


def validation_against_official(
    cohort: str,
    score_frame: pd.DataFrame,
    old_scores: pd.DataFrame,
    signature_name: str,
) -> pd.DataFrame:
    old = old_scores.set_index("cell_id").loc[score_frame.index]
    comparisons = [
        (f"ucell_{signature_name}", "UCell_VM_anchor"),
        (f"aucell_{signature_name}", "AUCell_VM_anchor"),
        (f"module_{signature_name}", "AddModuleScore_VM_anchor"),
    ]
    rows = []
    for new_field, old_field in comparisons:
        rho, p_value = stats.spearmanr(score_frame[new_field], old[old_field])
        rows.append(
            {
                "cohort": cohort,
                "python_field": new_field,
                "official_or_previous_field": old_field,
                "n_cells": len(score_frame),
                "spearman_rho": rho,
                "p_value": p_value,
                "max_absolute_difference": float(np.max(np.abs(score_frame[new_field] - old[old_field]))),
                "mean_absolute_difference": float(np.mean(np.abs(score_frame[new_field] - old[old_field]))),
            }
        )
    return pd.DataFrame(rows)


def run_cohort(cohort: str, config: dict[str, object]) -> dict[str, pd.DataFrame]:
    obj = ad.read_h5ad(config["h5ad"])
    restore_audit = pd.DataFrame()
    if cohort == "GSE138709":
        obj, restore_audit = restore_gse138709_anchor_columns(obj)
        restore_audit.to_csv(OUT / "GSE138709_restored_anchor_raw_UMI_audit.tsv", sep="\t", index=False)
    old_scores = pd.read_csv(config["old_scores"], sep="\t")
    selected = select_operational_cells(cohort, obj, old_scores)
    absent_anchors = [gene for gene in ANCHORS if gene not in selected.var_names]
    if absent_anchors:
        raise ValueError(f"{cohort}: anchors remain absent: {absent_anchors}")
    available_controls = [gene for gene in ENDOTHELIAL_CONTROLS if gene in selected.var_names]
    if len(available_controls) < 4:
        raise ValueError(f"{cohort}: too few measurable endothelial controls: {available_controls}")
    counts = selected.layers["counts"].tocsr()
    gene_index = pd.Series(np.arange(selected.n_vars), index=selected.var_names.astype(str))
    targets = ANCHORS + available_controls
    target_indices = [int(gene_index[gene]) for gene in targets]
    u_ranks, a_ranks = signature_ranks(counts, target_indices, seed=SEED + sum(map(ord, cohort)))
    target_position = {gene: targets.index(gene) for gene in targets}
    signatures = {"full6": ANCHORS, "endothelial": available_controls}
    signatures.update({f"drop_{gene}": [item for item in ANCHORS if item != gene] for gene in ANCHORS})
    score_frame = module_scores(selected, signatures)
    for name, genes in signatures.items():
        positions = [target_position[gene] for gene in genes]
        score_frame[f"ucell_{name}"] = ucell_score(
            u_ranks, positions, max_rank=min(UCELL_MAX_RANK, selected.n_vars)
        )
        score_frame[f"aucell_{name}"] = aucell_score(
            a_ranks, positions, auc_threshold=math.ceil(0.05 * selected.n_vars)
        )
    labels_80 = {name: call_high(score_frame, name, 0.80) for name in signatures if name != "endothelial"}
    labels_90 = {name: call_high(score_frame, name, 0.90) for name in signatures if name != "endothelial"}
    patient_field = str(config["patient"])
    patients = selected.obs[patient_field].astype(str).reset_index(drop=True)
    patients.index = score_frame.index
    classification_80 = classification_summary(cohort, labels_80, patients)
    classification_80["VM_score_quantile"] = 0.80
    classification_90 = classification_summary(cohort, labels_90, patients)
    classification_90["VM_score_quantile"] = 0.90
    threshold = pd.concat([classification_80, classification_90], ignore_index=True)
    corrected_full = labels_80["full6"]
    old_label = old_scores.set_index("cell_id").loc[score_frame.index, "VM_high_candidate"].astype(bool)
    old_vs_corrected = classification_summary(
        cohort,
        {"full6": corrected_full, "previous_reported_label": old_label},
        patients,
    )
    held_detail, held_summary = held_out_pseudobulk(cohort, selected, labels_80, patients)
    validation_signature = "full6" if cohort == "GSE125449" else "drop_MMP2"
    # For GSE138709 the previous official score used only CDH5, EPHA2 and MCAM.
    if cohort == "GSE138709":
        old_three = ["CDH5", "EPHA2", "MCAM"]
        signatures["old3"] = old_three
        positions = [target_position[gene] for gene in old_three]
        score_frame["ucell_old3"] = ucell_score(u_ranks, positions, min(UCELL_MAX_RANK, selected.n_vars))
        score_frame["aucell_old3"] = aucell_score(a_ranks, positions, math.ceil(0.05 * selected.n_vars))
        old3_module = module_scores(selected, {"old3": old_three})
        score_frame["module_old3"] = old3_module["module_old3"]
        validation_signature = "old3"
    validation = validation_against_official(cohort, score_frame, old_scores, validation_signature)
    export = score_frame.copy()
    export.insert(0, "cell_id", export.index)
    export.insert(1, "patient_id", patients.to_numpy())
    export["corrected_full6_high_q80"] = corrected_full.to_numpy()
    export["corrected_full6_high_q90"] = labels_90["full6"].to_numpy()
    export["previous_reported_high"] = old_label.to_numpy()
    for gene in ANCHORS:
        export[f"drop_{gene}_high_q80"] = labels_80[f"drop_{gene}"].to_numpy()
    export.to_csv(OUT / f"{cohort}_corrected_full6_and_leave_one_anchor_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    (OUT / f"{cohort}_measurable_endothelial_controls.txt").write_text(
        "\n".join(available_controls) + "\n", encoding="utf-8"
    )
    held_detail.to_csv(OUT / f"{cohort}_held_out_anchor_patient_pseudobulk.tsv", sep="\t", index=False)
    return {
        "threshold": threshold,
        "old_vs_corrected": old_vs_corrected,
        "held_summary": held_summary,
        "validation": validation,
    }


def main() -> None:
    outputs: dict[str, list[pd.DataFrame]] = {
        "threshold": [],
        "old_vs_corrected": [],
        "held_summary": [],
        "validation": [],
    }
    for cohort, config in COHORTS.items():
        cohort_outputs = run_cohort(cohort, config)
        for key in outputs:
            outputs[key].append(cohort_outputs[key])
    combined = {key: pd.concat(frames, ignore_index=True) for key, frames in outputs.items()}
    combined["threshold"].to_csv(
        OUT / "Table_R3_1_corrected_threshold_and_leave_one_anchor_reclassification.tsv", sep="\t", index=False
    )
    combined["old_vs_corrected"].to_csv(
        OUT / "Table_R3_2_previous_vs_corrected_full6_classification.tsv", sep="\t", index=False
    )
    combined["held_summary"].to_csv(
        OUT / "Table_R3_3_held_out_anchor_patient_pseudobulk_tests.tsv", sep="\t", index=False
    )
    combined["validation"].to_csv(
        OUT / "Table_R3_4_python_port_validation_against_official_scores.tsv", sep="\t", index=False
    )
    manifest = {
        "generated_for": "Reviewer 3 revision",
        "seed_AUCell_random_ties": SEED,
        "seed_scanpy_module_score": MODULE_SEED,
        "anchors": ANCHORS,
        "endothelial_controls": ENDOTHELIAL_CONTROLS,
        "UCell": "Python port of UCell 2.10.1 average-tie ranking and u_stat formula; maxRank=1500.",
        "AUCell": "Python port of AUCell 1.28.0 ranking/AUC formula; ties shuffled reproducibly; aucMaxRank=ceiling(0.05*n_genes).",
        "module_score": "scanpy.tl.score_genes on non-integrated library-normalized expression; ctrl_size=50, n_bins=25.",
        "classification": "All three six-gene implementation scores strictly above the cohort-specific upper-tail threshold and none of the three endothelial-control scores strictly above its 90th percentile.",
        "thresholds": [0.80, 0.90],
        "anchor_ablation": "Scores, percentile cutoffs and classifications were recomputed de novo after deleting each anchor.",
        "held_out_test": "The omitted anchor was not used in classification; its raw-count patient-by-group pseudobulk log1p CPM was compared with paired two-sided Wilcoxon tests and within-cohort BH correction across six anchors.",
        "GSE138709_correction": "MMP2, MMP9 and LAMC2 were restored from raw per-library UMI matrices before all six-gene calculations. The retained feature-reduced h5ad had omitted these genes.",
        "claim_boundary": "Upper-tail labels are operational summaries of a continuous high-score programme, not evidence for a discrete state or VM morphology.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("\nCorrected threshold and leave-one-anchor summary")
    print(combined["threshold"].to_string(index=False))
    print("\nPrevious versus corrected")
    print(combined["old_vs_corrected"].to_string(index=False))
    print("\nHeld-out anchor tests")
    print(combined["held_summary"].to_string(index=False))
    print("\nPython-port validation")
    print(combined["validation"].to_string(index=False))


if __name__ == "__main__":
    main()
