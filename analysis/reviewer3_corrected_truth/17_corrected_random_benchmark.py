from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(r"D:\第六篇大论文")
OUT = ROOT / "revision_round1_reviewer3_20261116/analysis_08_corrected_random_benchmark"
OUT.mkdir(parents=True, exist_ok=True)
ANCHORS = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
SEED = 20261116
N_RANDOM = 1000

COHORTS = {
    "GSE138709": {
        "h5ad": ROOT / "results/11_gse138709_strict_qc_map/gse138709_strict_qc_map.h5ad",
        "old": ROOT / "results/15_gse138709_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "corrected": ROOT / "revision_round1_reviewer3_20261116/analysis_01_corrected_scoring_and_anchor_ablation/GSE138709_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
    },
    "GSE125449": {
        "h5ad": ROOT / "results/43_gse125449_icca_malignant_epithelial_definition/gse125449_icca_operational_malignant_epithelium.h5ad",
        "old": ROOT / "results/46_gse125449_icca_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
        "corrected": ROOT / "revision_round1_reviewer3_20261116/analysis_01_corrected_scoring_and_anchor_ablation/GSE125449_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
    },
}


def restore_gse138709(obj):
    missing = [g for g in ANCHORS if g not in obj.var_names]
    if not missing:
        return obj
    raw_root = ROOT / "results/01_gse138709_import_v2"
    restored = {g: np.zeros(obj.n_obs, dtype=np.float32) for g in missing}
    pos = pd.Series(np.arange(obj.n_obs), index=obj.obs_names.astype(str))
    for d in sorted(x for x in raw_root.iterdir() if x.is_dir()):
        q, g, m = d / "cell_qc_unfiltered.tsv", d / "genes.tsv", d / "counts_cells_by_genes.npz"
        if not (q.exists() and g.exists() and m.exists()):
            continue
        qc = pd.read_csv(q, sep="\t")
        genes = pd.read_csv(g, sep="\t").iloc[:, 0].astype(str)
        matrix = sparse.load_npz(m).tocsr()
        common = qc.cell_id.astype(str).isin(pos.index)
        raw_rows = np.flatnonzero(common.to_numpy())
        target = pos.loc[qc.loc[common, "cell_id"].astype(str)].to_numpy(int)
        lookup = pd.Series(np.arange(len(genes)), index=genes)
        for gene in missing:
            if gene in lookup.index:
                restored[gene][target] = np.asarray(matrix[raw_rows, int(lookup[gene])].todense()).ravel()
    counts = obj.layers["counts"]
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    counts = sparse.hstack([counts, sparse.csr_matrix(np.column_stack([restored[g] for g in missing]))], format="csr")
    totals = np.asarray(counts.sum(axis=1)).ravel(); totals[totals <= 0] = 1
    x = sparse.diags(1e4/totals).dot(counts); x.data = np.log1p(x.data)
    new = ad.AnnData(X=x, obs=obj.obs.copy(), var=pd.DataFrame(index=list(obj.var_names.astype(str))+missing))
    new.layers["counts"] = counts
    return new


def paired_g(high, low):
    diff = np.asarray(high) - np.asarray(low)
    if len(diff) < 2 or np.std(diff, ddof=1) == 0:
        return np.nan
    correction = 1 - 3/(4*(len(diff)-1)-1) if len(diff) > 2 else 1
    return correction*np.mean(diff)/np.std(diff, ddof=1)


def run(cohort, cfg):
    rng = np.random.default_rng(SEED + sum(map(ord, cohort)))
    obj = ad.read_h5ad(cfg["h5ad"])
    if cohort == "GSE138709":
        obj = restore_gse138709(obj)
    old = pd.read_csv(cfg["old"], sep="\t")
    corrected = pd.read_csv(cfg["corrected"], sep="\t")
    obj = obj[old.cell_id.astype(str).to_numpy()].copy()
    labels = corrected.set_index("cell_id").loc[obj.obs_names, "corrected_full6_high_q80"].astype(bool)
    counts = obj.layers["counts"]
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    total = np.asarray(counts.sum(axis=1)).ravel(); total[total <= 0] = 1
    x = sparse.diags(1e4/total).dot(counts); x.data = np.log1p(x.data)
    mean = np.asarray(x.mean(axis=0)).ravel()
    detect = np.asarray((counts > 0).mean(axis=0)).ravel()
    metrics = pd.DataFrame({"gene": obj.var_names.astype(str), "mean": mean, "detect": detect}).set_index("gene")
    metrics["mean_rank"] = metrics["mean"].rank(pct=True)
    metrics["detect_rank"] = metrics["detect"].rank(pct=True)
    pools = {}
    for anchor in ANCHORS:
        a = metrics.loc[anchor]
        distance = np.sqrt((metrics.mean_rank-a.mean_rank)**2+(metrics.detect_rank-a.detect_rank)**2)
        pools[anchor] = distance.drop(index=[g for g in ANCHORS if g in distance.index]).sort_values().head(200).index.tolist()

    patients = obj.obs.patient_id.astype(str)
    patient_ids = []
    pb = []
    for patient in sorted(patients.unique()):
        for group, value in [("high", True), ("remainder", False)]:
            mask = patients.eq(patient).to_numpy() & labels.eq(value).to_numpy()
            if mask.sum() == 0:
                continue
            pb.append(np.asarray(counts[mask].sum(axis=0)).ravel())
            patient_ids.append((patient, group))
    pb = np.vstack(pb).T
    cpm = pb/np.maximum(pb.sum(axis=0, keepdims=True), 1)*1e6
    log = np.log1p(cpm)
    z = (log-log.mean(axis=1, keepdims=True))/np.where(log.std(axis=1, keepdims=True)==0, np.nan, log.std(axis=1, keepdims=True))
    frame = pd.DataFrame(z, index=obj.var_names.astype(str), columns=pd.MultiIndex.from_tuples(patient_ids))
    eligible = [p for p in sorted(patients.unique()) if (p, "high") in frame.columns and (p, "remainder") in frame.columns]

    def effect(genes):
        score = frame.loc[genes].mean(axis=0)
        return paired_g([score[(p, "high")] for p in eligible], [score[(p, "remainder")] for p in eligible])

    actual = effect(ANCHORS)
    null = []
    diagnostics = []
    for iteration in range(1, N_RANDOM+1):
        chosen = []
        for anchor in ANCHORS:
            available = [g for g in pools[anchor][:100] if g not in chosen]
            gene = rng.choice(available)
            chosen.append(gene)
            diagnostics.append({"cohort": cohort, "iteration": iteration, "anchor": anchor, "matched_gene": gene})
        null.append({"cohort": cohort, "iteration": iteration, "paired_Hedges_g": effect(chosen), "genes": ";".join(chosen)})
    null = pd.DataFrame(null)
    finite = null.paired_Hedges_g.dropna()
    empirical = (1+(finite.abs() >= abs(actual)).sum())/(1+len(finite))
    summary = {"cohort": cohort, "n_paired_patients": len(eligible), "n_random_sets": len(finite),
               "six_anchor_paired_Hedges_g": actual, "null_mean": finite.mean(), "null_sd": finite.std(ddof=1),
               "empirical_two_sided_p": empirical,
               "interpretation": "Internal calibration only because the corrected labels are anchor-defined."}
    return summary, null, pd.DataFrame(diagnostics)


def main():
    summaries, nulls, diagnostics = [], [], []
    for cohort, cfg in COHORTS.items():
        s, n, d = run(cohort, cfg)
        summaries.append(s); nulls.append(n); diagnostics.append(d)
    pd.DataFrame(summaries).to_csv(OUT / "Table_U12_corrected_detection_expression_random_benchmark.tsv", sep="\t", index=False)
    pd.concat(nulls).to_csv(OUT / "corrected_random_effect_distribution.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.concat(diagnostics).to_csv(OUT / "corrected_matching_diagnostics.tsv.gz", sep="\t", index=False, compression="gzip")
    (OUT / "analysis_manifest.json").write_text(json.dumps({"seed": SEED, "n_random_sets": N_RANDOM,
        "source_of_truth": "Reviewer 3 corrected complete-six-gene q80 labels",
        "boundary": "Circular internal calibration, not specificity or independent validation."}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
