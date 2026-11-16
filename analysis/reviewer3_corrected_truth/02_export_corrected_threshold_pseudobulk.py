from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(r"D:\第六篇大论文")
REVISION = ROOT / "revision_round1_reviewer3_20261116"
SCORING = REVISION / "analysis_01_corrected_scoring_and_anchor_ablation"
OUT = REVISION / "analysis_02_corrected_threshold_pseudobulk"
INPUTS = OUT / "pseudobulk_inputs"
OUT.mkdir(parents=True, exist_ok=True)
INPUTS.mkdir(parents=True, exist_ok=True)
MIN_HIGH_CELLS = 5


def load_scoring_module():
    script = REVISION / "01_recompute_six_gene_and_leave_one_anchor_labels.py"
    spec = importlib.util.spec_from_file_location("reviewer3_scoring", script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_gse138709() -> tuple[ad.AnnData, pd.DataFrame, str, str]:
    module = load_scoring_module()
    obj = ad.read_h5ad(ROOT / "results/11_gse138709_strict_qc_map/gse138709_strict_qc_map.h5ad")
    obj, _ = module.restore_gse138709_anchor_columns(obj)
    scores = pd.read_csv(
        SCORING / "GSE138709_corrected_full6_and_leave_one_anchor_scores.tsv.gz", sep="\t"
    ).set_index("cell_id")
    selected = obj[scores.index].copy()
    subclusters = ad.read_h5ad(
        ROOT / "results/30_gse138709_malignant_epithelial_subclusters/gse138709_malignant_epithelial_subclusters.h5ad",
        backed="r",
    )
    subcluster_values = subclusters.obs["malignant_epithelial_subcluster"].astype(str)
    subclusters.file.close()
    selected.obs["malignant_epithelial_subcluster"] = subcluster_values.loc[selected.obs_names].to_numpy()
    core = selected.obs["malignant_epithelial_subcluster"].astype(str).isin({"0", "1", "2", "3"})
    selected = selected[core.to_numpy()].copy()
    scores = scores.loc[selected.obs_names]
    return selected, scores, "sample_id", "patient_id"


def prepare_gse125449() -> tuple[ad.AnnData, pd.DataFrame, str, str]:
    obj = ad.read_h5ad(
        ROOT / "results/43_gse125449_icca_malignant_epithelial_definition/gse125449_icca_operational_malignant_epithelium.h5ad"
    )
    scores = pd.read_csv(
        SCORING / "GSE125449_corrected_full6_and_leave_one_anchor_scores.tsv.gz", sep="\t"
    ).set_index("cell_id")
    selected = obj[scores.index].copy()
    scores = scores.loc[selected.obs_names]
    return selected, scores, "library_id", "patient_id"


def export_one(
    cohort: str,
    obj: ad.AnnData,
    scores: pd.DataFrame,
    library_field: str,
    patient_field: str,
    quantile: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    label_field = f"corrected_full6_high_q{quantile}"
    labels = scores[label_field].astype(bool)
    metadata = obj.obs[[library_field, patient_field]].copy()
    metadata["high_score_group"] = labels.loc[obj.obs_names].to_numpy()
    metadata["score_group"] = np.where(metadata["high_score_group"], "high", "low")
    eligibility = (
        metadata.groupby([patient_field, library_field, "score_group"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .pivot_table(index=[patient_field, library_field], columns="score_group", values="n_cells", fill_value=0)
        .reset_index()
        .rename_axis(columns=None)
    )
    for group in ["high", "low"]:
        if group not in eligibility:
            eligibility[group] = 0
    eligibility["eligible"] = eligibility["high"].ge(MIN_HIGH_CELLS) & eligibility["low"].gt(0)
    eligibility["exclusion_reason"] = np.where(
        eligibility["eligible"],
        "",
        np.where(eligibility["high"].eq(0), "no_high_score_cells", f"high_score_cells_below_{MIN_HIGH_CELLS}"),
    )
    included = set(eligibility.loc[eligibility["eligible"], library_field].astype(str))
    selected_mask = metadata[library_field].astype(str).isin(included)
    selected = obj[selected_mask.to_numpy()].copy()
    selected_meta = metadata.loc[selected.obs_names]
    counts = selected.layers["counts"].tocsr()
    vectors = []
    rows = []
    for (patient, library, group), frame in selected_meta.groupby(
        [patient_field, library_field, "score_group"], observed=True, sort=True
    ):
        positions = selected.obs_names.get_indexer(frame.index)
        vectors.append(sparse.csr_matrix(counts[positions].sum(axis=0)))
        rows.append(
            {
                "patient_id": str(patient),
                "library_id": str(library),
                "score_group": group,
                "n_cells": len(positions),
                "pseudobulk_id": f"{library}__{group}",
            }
        )
    sample_table = pd.DataFrame(rows)
    paired = sample_table.groupby("library_id", observed=True)["score_group"].nunique()
    if len(sample_table) and not paired.eq(2).all():
        raise ValueError(f"{cohort} q{quantile}: an included library is not paired")
    output_dir = INPUTS / f"{cohort}_q{quantile}"
    output_dir.mkdir(parents=True, exist_ok=True)
    if vectors:
        matrix = sparse.vstack(vectors).transpose().tocsr()
        with (output_dir / "pseudobulk_counts_genes_by_samples.mtx").open("wb") as handle:
            mmwrite(handle, sparse.csc_matrix(matrix), field="integer")
    pd.DataFrame({"gene_symbol": selected.var_names.astype(str)}).to_csv(
        output_dir / "genes.tsv", sep="\t", index=False
    )
    sample_table.to_csv(output_dir / "pseudobulk_samples.tsv", sep="\t", index=False)
    eligibility.to_csv(output_dir / "library_eligibility.tsv", sep="\t", index=False)
    summary = {
        "cohort": cohort,
        "VM_score_quantile": quantile / 100,
        "minimum_high_cells": MIN_HIGH_CELLS,
        "eligible_libraries": int(sample_table["library_id"].nunique()) if len(sample_table) else 0,
        "eligible_patients": int(sample_table["patient_id"].nunique()) if len(sample_table) else 0,
        "pseudobulks": len(sample_table),
        "minimum_high_cells_among_eligible": int(
            sample_table.loc[sample_table["score_group"].eq("high"), "n_cells"].min()
        ) if len(sample_table) else 0,
        "input_directory": str(output_dir.relative_to(ROOT)),
    }
    return eligibility.assign(cohort=cohort, VM_score_quantile=quantile / 100), summary


def main() -> None:
    cohort_objects = {
        "GSE138709": prepare_gse138709(),
        "GSE125449": prepare_gse125449(),
    }
    eligibility_tables = []
    summaries = []
    for cohort, (obj, scores, library_field, patient_field) in cohort_objects.items():
        for quantile in [80, 90]:
            eligibility, summary = export_one(
                cohort, obj, scores, library_field, patient_field, quantile
            )
            eligibility_tables.append(eligibility)
            summaries.append(summary)
    pd.concat(eligibility_tables, ignore_index=True).to_csv(
        OUT / "Table_R3_5_corrected_threshold_pseudobulk_eligibility.tsv", sep="\t", index=False
    )
    pd.DataFrame(summaries).to_csv(
        OUT / "corrected_threshold_pseudobulk_summary.tsv", sep="\t", index=False
    )
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "labels": "Corrected complete-six-gene classifications from Reviewer 3 analysis 01.",
                "thresholds": [0.80, 0.90],
                "minimum_high_cells_per_library": MIN_HIGH_CELLS,
                "GSE138709_scope": "Operational malignant epithelial subclusters 0-3, matching the original core-epithelial SecAct input.",
                "GSE125449_scope": "Operational malignant epithelial iCCA cells only.",
                "unit": "Paired library-level high-score and low-score raw-count pseudobulk.",
                "boundary": "These dichotomies are sensitivity summaries of a continuous score and are not discrete biological states.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
