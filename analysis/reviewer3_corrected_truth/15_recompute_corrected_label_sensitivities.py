from __future__ import annotations

import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


ROOT = Path(r"D:\第六篇大论文")
R3 = ROOT / "revision_round1_reviewer3_20261116"
OUT = R3 / "analysis_06_unified_corrected_sensitivities"
OUT.mkdir(parents=True, exist_ok=True)

ANCHORS = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
ENDOTHELIAL_CONTROLS = ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "EGFL7", "PLVAP", "CA4"]
HIGH_QUANTILES = [0.70, 0.75, 0.80, 0.85, 0.90]
ENDOTHELIAL_QUANTILES = [0.80, 0.85, 0.90, 0.95]
MIN_HIGH_CELLS = 5
MODULE_SEED = 20260830
UCELL_MAX_RANK = 1500

COHORTS = {
    "GSE138709": {
        "h5ad": ROOT / "results/11_gse138709_strict_qc_map/gse138709_strict_qc_map.h5ad",
        "old_scores": ROOT / "results/15_gse138709_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "corrected_scores": R3 / "analysis_01_corrected_scoring_and_anchor_ablation/GSE138709_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
        "replicate": "sample_id",
    },
    "GSE125449": {
        "h5ad": ROOT / "results/43_gse125449_icca_malignant_epithelial_definition/gse125449_icca_operational_malignant_epithelium.h5ad",
        "old_scores": ROOT / "results/46_gse125449_icca_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "corrected_scores": R3 / "analysis_01_corrected_scoring_and_anchor_ablation/GSE125449_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
        "replicate": "library_id",
    },
}


def restore_gse138709(obj: ad.AnnData) -> ad.AnnData:
    missing = [g for g in ANCHORS if g not in obj.var_names]
    if not missing:
        return obj
    raw_root = ROOT / "results/01_gse138709_import_v2"
    restored = {g: np.zeros(obj.n_obs, dtype=np.float32) for g in missing}
    positions = pd.Series(np.arange(obj.n_obs), index=obj.obs_names.astype(str))
    for sample_dir in sorted(p for p in raw_root.iterdir() if p.is_dir()):
        qc_path = sample_dir / "cell_qc_unfiltered.tsv"
        genes_path = sample_dir / "genes.tsv"
        matrix_path = sample_dir / "counts_cells_by_genes.npz"
        if not (qc_path.exists() and genes_path.exists() and matrix_path.exists()):
            continue
        qc = pd.read_csv(qc_path, sep="\t")
        genes = pd.read_csv(genes_path, sep="\t").iloc[:, 0].astype(str)
        matrix = sparse.load_npz(matrix_path).tocsr()
        common = qc["cell_id"].astype(str).isin(positions.index)
        if not common.any():
            continue
        raw_rows = np.flatnonzero(common.to_numpy())
        target_rows = positions.loc[qc.loc[common, "cell_id"].astype(str)].to_numpy(int)
        gene_index = pd.Series(np.arange(len(genes)), index=genes)
        for gene in missing:
            if gene in gene_index.index:
                restored[gene][target_rows] = np.asarray(
                    matrix[raw_rows, int(gene_index[gene])].todense()
                ).ravel()
    counts = obj.layers["counts"]
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    all_counts = sparse.hstack(
        [counts, sparse.csr_matrix(np.column_stack([restored[g] for g in missing]))],
        format="csr",
    )
    totals = np.asarray(all_counts.sum(axis=1)).ravel()
    totals[totals <= 0] = 1
    normalized = sparse.diags(1e4 / totals).dot(all_counts)
    normalized.data = np.log1p(normalized.data)
    out = ad.AnnData(
        X=normalized.tocsr(),
        obs=obj.obs.copy(),
        var=pd.DataFrame(index=list(obj.var_names.astype(str)) + missing),
    )
    out.layers["counts"] = all_counts
    return out


def select_operational(obj: ad.AnnData, old_scores: pd.DataFrame) -> ad.AnnData:
    ids = old_scores["cell_id"].astype(str).to_numpy()
    if not pd.Index(ids).isin(obj.obs_names.astype(str)).all():
        raise ValueError("Operational cell identities do not match source object")
    return obj[ids].copy()


def call_from_export(scores: pd.DataFrame, high_q: float, endo_q: float) -> tuple[pd.Series, pd.Series, pd.Series]:
    high = pd.concat(
        [scores[c].gt(scores[c].quantile(high_q)) for c in ["module_full6", "ucell_full6", "aucell_full6"]],
        axis=1,
    ).all(axis=1)
    endo = pd.concat(
        [scores[c].gt(scores[c].quantile(endo_q)) for c in ["module_endothelial", "ucell_endothelial", "aucell_endothelial"]],
        axis=1,
    ).any(axis=1)
    return high, endo, high & ~endo


def jaccard(a: pd.Series, b: pd.Series) -> float:
    av, bv = a.to_numpy(bool), b.to_numpy(bool)
    union = np.logical_or(av, bv).sum()
    return float(np.logical_and(av, bv).sum() / union) if union else 1.0


def kappa(a: pd.Series, b: pd.Series) -> float:
    av, bv = a.to_numpy(bool), b.to_numpy(bool)
    obs = np.mean(av == bv)
    exp = np.mean(av) * np.mean(bv) + (1 - np.mean(av)) * (1 - np.mean(bv))
    return float((obs - exp) / (1 - exp)) if exp < 1 else np.nan


def parse_enrichr(path: Path) -> dict[str, list[str]]:
    out = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                out[fields[0]] = sorted({g for g in fields[2:] if g})
    return out


def comparator_sets() -> dict[str, list[str]]:
    root = ROOT / "revision_round1_reviewer1_20261115/reference_gene_sets"
    hallmark = parse_enrichr(root / "MSigDB_Hallmark_2020.txt")
    reactome = parse_enrichr(root / "Reactome_2022.txt")
    go = parse_enrichr(root / "GO_Biological_Process_2025.txt")
    panglao = parse_enrichr(root / "PanglaoDB_Augmented_2021.txt")
    return {
        "EMT_Hallmark_anchor_excluded": sorted(set(hallmark["Epithelial Mesenchymal Transition"]) - set(ANCHORS)),
        "Hypoxia_Hallmark_anchor_excluded": sorted(set(hallmark["Hypoxia"]) - set(ANCHORS)),
        "ECM_Reactome_anchor_excluded": sorted(set(reactome["Extracellular Matrix Organization R-HSA-1474244"]) - set(ANCHORS)),
        "Epithelial_migration_GO_anchor_excluded": sorted(set(go["Positive Regulation of Epithelial Cell Migration (GO:0010634)"]) - set(ANCHORS)),
        "Endothelial_identity_Panglao_independent": sorted(
            set(panglao["Endothelial Cells"]) - set(ANCHORS) - {"PECAM1", "VWF", "KDR", "EMCN"}
        ),
    }


def mean_score(matrix, var_names, genes):
    lookup = {str(g): i for i, g in enumerate(var_names)}
    present = [g for g in genes if g in lookup]
    if not present:
        return np.full(matrix.shape[0], np.nan), []
    values = matrix[:, [lookup[g] for g in present]].mean(axis=1)
    return np.asarray(values).ravel(), present


def hedges_g(values, labels):
    high, low = values[labels], values[~labels]
    high, low = high[np.isfinite(high)], low[np.isfinite(low)]
    if len(high) < 2 or len(low) < 2:
        return np.nan, np.nan, np.nan
    pooled = math.sqrt(((len(high)-1)*np.var(high, ddof=1)+(len(low)-1)*np.var(low, ddof=1))/(len(high)+len(low)-2))
    if pooled == 0:
        return np.nan, np.nan, np.nan
    raw = (np.mean(high)-np.mean(low))/pooled
    g = (1 - 3/(4*(len(high)+len(low))-9))*raw
    var = (len(high)+len(low))/(len(high)*len(low)) + g*g/(2*(len(high)+len(low)-2))
    margin = stats.t.ppf(.975, len(high)+len(low)-2)*math.sqrt(var)
    return g, g-margin, g+margin


def main():
    threshold_rows, replicate_rows, loo_rows = [], [], []
    programme_rows, patient_programme_rows, endothelial_rows = [], [], []
    score_corr_rows, qc_summary_rows, qc_test_rows, mito_rows = [], [], [], []
    programme_definitions = comparator_sets()

    for cohort, cfg in COHORTS.items():
        old_scores = pd.read_csv(cfg["old_scores"], sep="\t")
        scores = pd.read_csv(cfg["corrected_scores"], sep="\t")
        if not old_scores["cell_id"].astype(str).equals(scores["cell_id"].astype(str)):
            raise ValueError(f"{cohort}: corrected scores are not in the frozen operational-cell order")
        scores[cfg["replicate"]] = old_scores[cfg["replicate"]].astype(str).to_numpy()
        primary = scores["corrected_full6_high_q80"].astype(bool)
        _, _, reproduced = call_from_export(scores, .80, .90)
        if not reproduced.equals(primary):
            raise ValueError(f"{cohort}: corrected q80 export does not reproduce")

        for hq in HIGH_QUANTILES:
            for eq in ENDOTHELIAL_QUANTILES:
                triple, endo, high = call_from_export(scores, hq, eq)
                frame = scores[["patient_id", cfg["replicate"]]].copy()
                frame["high"] = high
                frame["low"] = ~high
                by_rep = frame.groupby(["patient_id", cfg["replicate"]], observed=True).agg(
                    high_cells=("high", "sum"), low_cells=("low", "sum")
                ).reset_index()
                eligible = by_rep.high_cells.ge(MIN_HIGH_CELLS) & by_rep.low_cells.gt(0)
                by_rep.insert(0, "cohort", cohort)
                by_rep.insert(1, "high_quantile", hq)
                by_rep.insert(2, "endothelial_quantile", eq)
                by_rep["eligible_for_paired_pseudobulk"] = eligible
                replicate_rows.append(by_rep)
                threshold_rows.append({
                    "cohort": cohort, "high_quantile": hq, "endothelial_quantile": eq,
                    "malignant_epithelial_cells": len(scores), "triple_consensus_cells": int(triple.sum()),
                    "endothelial_flagged_cells": int(endo.sum()),
                    "endothelial_flagged_triple_cells": int((triple & endo).sum()),
                    "vm_high_cells": int(high.sum()), "vm_high_fraction": float(high.mean()),
                    "patients_total": int(scores.patient_id.nunique()),
                    "patients_with_at_least_1_vm_high": int(scores.loc[high, "patient_id"].nunique()),
                    "eligible_replicates": int(eligible.sum()),
                    "eligible_patients": int(by_rep.loc[eligible, "patient_id"].nunique()),
                    "jaccard_with_corrected_q80": jaccard(high, primary),
                    "kappa_with_corrected_q80": kappa(high, primary),
                })

        for omitted in sorted(scores.patient_id.astype(str).unique()):
            retained = scores.patient_id.astype(str).ne(omitted)
            sub = scores.loc[retained].reset_index(drop=True)
            _, endo, high = call_from_export(sub, .80, .90)
            ref = primary.loc[retained].reset_index(drop=True)
            frame = sub[["patient_id", cfg["replicate"]]].copy()
            frame["high"] = high
            frame["low"] = ~high
            by_rep = frame.groupby(["patient_id", cfg["replicate"]], observed=True).agg(
                high_cells=("high", "sum"), low_cells=("low", "sum")
            ).reset_index()
            eligible = by_rep.high_cells.ge(5) & by_rep.low_cells.gt(0)
            loo_rows.append({
                "cohort": cohort, "omitted_patient": omitted,
                "retained_patients": int(frame.patient_id.nunique()), "retained_cells": len(frame),
                "vm_high_cells": int(high.sum()), "vm_high_fraction": float(high.mean()),
                "endothelial_flagged_cells": int(endo.sum()),
                "eligible_replicates": int(eligible.sum()),
                "eligible_patients": int(by_rep.loc[eligible, "patient_id"].nunique()),
                "jaccard_with_corrected_full_cohort_q80_on_retained_cells": jaccard(high, ref),
                "kappa_with_corrected_full_cohort_q80_on_retained_cells": kappa(high, ref),
            })

        for a, b in [("module_full6", "ucell_full6"), ("module_full6", "aucell_full6"), ("ucell_full6", "aucell_full6")]:
            rho, p = stats.spearmanr(scores[a], scores[b])
            score_corr_rows.append({"cohort": cohort, "score_1": a, "score_2": b,
                                    "n_cells": len(scores), "spearman_rho": rho, "p_value": p})

        obj = ad.read_h5ad(cfg["h5ad"])
        if cohort == "GSE138709":
            obj = restore_gse138709(obj)
        obj = select_operational(obj, old_scores)
        counts = obj.layers["counts"]
        counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
        matrix = obj.X.tocsr() if sparse.issparse(obj.X) else sparse.csr_matrix(obj.X)
        metadata = obj.obs.copy()
        metadata["patient_id"] = scores.patient_id.astype(str).to_numpy()
        metadata["high"] = primary.to_numpy()

        q = metadata[["patient_id", "n_counts", "n_genes", "pct_mito", "high"]].copy()
        q["score_group"] = np.where(q.high, "high", "remainder")
        for (pid, group), sub in q.groupby(["patient_id", "score_group"], observed=True):
            for metric in ["n_counts", "n_genes", "pct_mito"]:
                vals = sub[metric].astype(float)
                qc_summary_rows.append({"cohort": cohort, "patient_id": pid, "score_group": group,
                                        "metric": metric, "n_cells": len(vals), "median": vals.median(),
                                        "mean": vals.mean(), "q25": vals.quantile(.25), "q75": vals.quantile(.75)})
        summary = pd.DataFrame([x for x in qc_summary_rows if x["cohort"] == cohort])
        for metric in ["n_counts", "n_genes", "pct_mito"]:
            wide = summary[summary.metric.eq(metric)].pivot(index="patient_id", columns="score_group", values="median").dropna()
            if len(wide) >= 2 and (wide.high - wide.remainder).abs().sum() > 0:
                w, p = stats.wilcoxon(wide.high, wide.remainder)
            else:
                w, p = np.nan, np.nan
            qc_test_rows.append({"cohort": cohort, "metric": metric, "n_paired_patients": len(wide),
                                 "paired_patient_Wilcoxon_W": w, "paired_patient_p": p,
                                 "median_patient_difference_high_minus_remainder": float((wide.high-wide.remainder).median()) if len(wide) else np.nan})
        for ceiling in [15.0, 20.0]:
            kept = q.pct_mito.le(ceiling)
            for group in ["high", "remainder"]:
                base = q.score_group.eq(group)
                mito_rows.append({"cohort": cohort, "mitochondrial_ceiling_percent": ceiling,
                                  "score_group": group, "n_original": int(base.sum()),
                                  "n_retained": int((base & kept).sum()),
                                  "retained_fraction": float((base & kept).sum()/base.sum()) if base.sum() else np.nan})
            mito_rows.append({"cohort": cohort, "mitochondrial_ceiling_percent": ceiling,
                              "score_group": "all_malignant", "n_original": len(q),
                              "n_retained": int(kept.sum()), "retained_fraction": float(kept.mean()),
                              "high_score_fraction_after_filter": float(q.loc[kept, "high"].mean())})

        anchor_score, _ = mean_score(matrix, obj.var_names, ANCHORS)
        programme_frame = pd.DataFrame(index=scores.index)
        for name, genes in programme_definitions.items():
            values, present = mean_score(matrix, obj.var_names, genes)
            programme_frame[name] = values
            g, low, high_ci = hedges_g(values, primary.to_numpy())
            programme_rows.append({
                "cohort": cohort, "program": name, "cells": len(values),
                "patients": int(scores.patient_id.nunique()), "high_score_cells": int(primary.sum()),
                "remainder_cells": int((~primary).sum()), "genes_present": len(present),
                "spearman_with_six_anchor_expression_score": stats.spearmanr(values, anchor_score, nan_policy="omit").statistic,
                "hedges_g_high_vs_remainder": g, "hedges_g_ci95_low": low, "hedges_g_ci95_high": high_ci,
            })
            for patient in sorted(scores.patient_id.astype(str).unique()):
                mask = scores.patient_id.astype(str).eq(patient).to_numpy()
                labels = primary.to_numpy()[mask]
                patient_values = values[mask]
                pg, plow, phigh = hedges_g(patient_values, labels)
                patient_programme_rows.append({
                    "cohort": cohort, "patient_id": patient, "program": name,
                    "high_score_cells": int(labels.sum()), "remainder_cells": int((~labels).sum()),
                    "mean_difference_high_minus_remainder": float(np.nanmean(patient_values[labels])-np.nanmean(patient_values[~labels])) if labels.any() and (~labels).any() else np.nan,
                    "hedges_g": pg, "ci95_low": plow, "ci95_high": phigh,
                })
        endothelial = programme_frame["Endothelial_identity_Panglao_independent"]
        for quantile in [.80, .85, .90, .95]:
            cutoff = float(endothelial.quantile(quantile))
            retained = primary & endothelial.le(cutoff)
            endothelial_rows.append({
                "cohort": cohort, "independent_endothelial_exclusion_quantile": quantile,
                "independent_endothelial_score_cutoff": cutoff,
                "original_high_score_cells": int(primary.sum()),
                "high_score_cells_after_independent_endothelial_exclusion": int(retained.sum()),
                "fraction_original_high_score_retained": float(retained.sum()/primary.sum()),
                "patients_with_retained_high_score_cells": int(scores.loc[retained, "patient_id"].nunique()),
            })

    pd.DataFrame(threshold_rows).to_csv(OUT / "Table_U1_corrected_threshold_grid.tsv", sep="\t", index=False)
    pd.concat(replicate_rows, ignore_index=True).to_csv(OUT / "Table_U2_corrected_threshold_replicate_counts.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT / "Table_U3_corrected_leave_one_patient_out.tsv", sep="\t", index=False)
    pd.DataFrame(score_corr_rows).to_csv(OUT / "Table_U4_corrected_score_correlations.tsv", sep="\t", index=False)
    pd.DataFrame(qc_summary_rows).to_csv(OUT / "Table_U5_corrected_patient_QC_summaries.tsv", sep="\t", index=False)
    pd.DataFrame(qc_test_rows).to_csv(OUT / "Table_U6_corrected_patient_QC_tests.tsv", sep="\t", index=False)
    pd.DataFrame(mito_rows).to_csv(OUT / "Table_U7_corrected_mitochondrial_ceiling_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(programme_rows).to_csv(OUT / "Table_U8_corrected_comparator_programmes.tsv", sep="\t", index=False)
    pd.DataFrame(patient_programme_rows).to_csv(OUT / "Table_U9_corrected_patient_programme_effects.tsv", sep="\t", index=False)
    pd.DataFrame(endothelial_rows).to_csv(OUT / "Table_U10_corrected_independent_endothelial_exclusion.tsv", sep="\t", index=False)
    definition_rows = [
        {"program": program, "gene": gene}
        for program, genes in programme_definitions.items()
        for gene in genes
    ]
    pd.DataFrame(definition_rows).to_csv(
        OUT / "Table_U0_comparator_programme_gene_definitions.tsv", sep="\t", index=False
    )
    (OUT / "analysis_manifest.json").write_text(json.dumps({
        "source_of_truth": "Reviewer 3 complete-six-gene corrected score exports",
        "high_threshold_grid": HIGH_QUANTILES,
        "endothelial_threshold_grid": ENDOTHELIAL_QUANTILES,
        "primary": "q80 complete-six-gene implementation consensus with q90 endothelial-control exclusion",
        "claim_boundary": "All binary analyses are operational upper-tail sensitivity audits of a continuous VM-related programme."
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
