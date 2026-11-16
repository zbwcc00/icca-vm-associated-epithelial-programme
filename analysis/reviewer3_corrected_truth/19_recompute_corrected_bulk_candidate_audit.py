from __future__ import annotations

import json
import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:\第六篇大论文")
R3 = ROOT / "revision_round1_reviewer3_20261116"
OUT = R3 / "analysis_10_corrected_bulk_candidate_audit"
OUT.mkdir(parents=True, exist_ok=True)

VM_ANCHORS = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
PAN_MYELOID = ["AIF1", "CSF1R", "FCER1G", "LST1", "PTPRC", "TYROBP", "LILRB1", "LILRB2"]
FIBROBLAST = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM", "PDGFRA"]
ENDOTHELIAL = ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "EGFL7", "PLVAP", "CA4"]


def standardize(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def score(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if f"member::{g}" in df.columns]
    if not present:
        return pd.Series(np.nan, index=df.index), []
    values = df[[f"member::{g}" for g in present]].mean(axis=1)
    return standardize(values), present


def score_matrix(matrix: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns), []
    return standardize(matrix.loc[present].mean(axis=0)), present


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=0).replace(0, np.nan), axis=0)


def read_gse107943() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = ROOT / "revision_round1_reviewer2_20261115/external_inputs/bulk/GSE107943_RPKM.txt.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = next(csv.reader([handle.readline().rstrip("\n")], delimiter="\t"))
        labels = next(csv.reader([handle.readline().rstrip("\n")], delimiter="\t"))
        expr = pd.read_csv(handle, sep="\t", names=header, header=None, dtype={"Genesymbol": str}, low_memory=False)
    mapping = {lab.strip(): col for lab, col in zip(labels[7:], header[7:]) if lab.strip().endswith("T")}
    meta = pd.read_csv(ROOT / "results/81_bulk_clinical_cohort_audit/GSE107943_sample_metadata_audited.tsv", sep="\t", dtype=str)
    meta = meta.loc[meta.audited_sample_role.eq("iCCA_tumor")].copy()
    meta["expression_column"] = meta.title.map(mapping)
    cols = meta.expression_column.tolist()
    values = expr[cols].apply(pd.to_numeric, errors="coerce")
    table = pd.concat([expr[["Genesymbol"]], values], axis=1)
    table["Genesymbol"] = table.Genesymbol.astype(str).str.upper().str.strip()
    table = table.loc[table.Genesymbol.str.fullmatch(r"[A-Z0-9-]+", na=False)].copy()
    table["variance"] = table[cols].var(axis=1, skipna=True)
    matrix = table.sort_values(["Genesymbol", "variance"], ascending=[True, False]).drop_duplicates("Genesymbol").set_index("Genesymbol")[cols]
    matrix = np.log2(matrix.clip(lower=0) + 0.1)
    matrix = zscore_rows(matrix)
    matrix.columns = meta.gsm_accession
    return matrix, meta.set_index("gsm_accession")


def read_gse26566() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_path = ROOT / "revision_round1_reviewer2_20261115/external_inputs/bulk/GSE26566_series_matrix.txt.gz"
    annot_path = R3 / "external_inputs/GPL6104.annot.gz"
    if not annot_path.exists():
        raise FileNotFoundError(f"Missing {annot_path}; download is handled before this script is rerun")
    meta = pd.read_csv(ROOT / "results/85_external_bulk_cohort_intake_audit/GSE26566_sample_metadata_audited.tsv", sep="\t", dtype=str, index_col="gsm_accession")
    samples = meta.index[meta.audited_sample_role.eq("cholangiocarcinoma_tumor")].tolist()
    with gzip.open(expr_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
        expr = pd.read_csv(handle, sep="\t", dtype={"ID_REF": str}, low_memory=False)
    expr = expr.loc[~expr.ID_REF.astype(str).str.startswith("!")].rename(columns={"ID_REF": "probe_id"}).set_index("probe_id")[samples].apply(pd.to_numeric, errors="coerce")
    expr = np.log2(expr.clip(lower=1.0))
    with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                break
        annot = pd.read_csv(handle, sep="\t", dtype=str, low_memory=False)
    annot = annot[["ID", "Gene symbol"]].rename(columns={"ID": "probe_id", "Gene symbol": "gene_symbol"})
    annot.gene_symbol = annot.gene_symbol.astype(str).str.upper().str.strip()
    annot.gene_symbol = annot.gene_symbol.str.split(r"\s*(?:///|;|,)\s*")
    annot = annot.explode("gene_symbol")
    annot = annot.loc[annot.gene_symbol.str.fullmatch(r"[A-Z0-9-]+", na=False) & annot.probe_id.isin(expr.index)].drop_duplicates()
    annot["variance"] = annot.probe_id.map(expr.var(axis=1, ddof=0))
    chosen = annot.sort_values(["gene_symbol", "variance", "probe_id"], ascending=[True, False, True]).drop_duplicates("gene_symbol")
    matrix = expr.loc[chosen.probe_id].copy()
    matrix.index = chosen.gene_symbol.to_numpy()
    matrix = zscore_rows(matrix)
    return matrix, meta.loc[samples]


def fit(df: pd.DataFrame, outcome: str, vm: str) -> dict[str, float | int]:
    covars = [vm, "pan_myeloid_total_proxy", "fibroblast_proxy", "endothelial_proxy"]
    d = df[[outcome, *covars]].dropna()
    if len(d) < 25:
        return {"n": len(d), "beta": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
    model = sm.OLS(d[outcome], sm.add_constant(d[covars])).fit(cov_type="HC3")
    ci = model.conf_int().loc[vm]
    return {"n": len(d), "beta": float(model.params[vm]), "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]), "p_value": float(model.pvalues[vm])}


def add_fdr(df: pd.DataFrame, p_col: str, out_col: str) -> None:
    ok = df[p_col].notna()
    df[out_col] = np.nan
    if ok.any():
        df.loc[ok, out_col] = multipletests(df.loc[ok, p_col], method="fdr_bh")[1]


constraint = pd.read_csv(R3 / "analysis_04_corrected_secact_integration/Table_R3_13_corrected_q80_source_and_external_constraint.tsv", sep="\t")
for col in ["nominal_positive_both", "external_RNA_compatible"]:
    constraint[col] = constraint[col].astype(str).str.lower().eq("true")
candidate = constraint.loc[constraint.nominal_positive_both & constraint.external_RNA_compatible].copy()
assert candidate.signature_group.nunique() == 26

member = pd.read_csv(ROOT / "results/78_secact_multilayer_evidence_matrix/secact_source_member_multilayer_evidence_matrix.tsv", sep="\t")
member["external_state_restricted_source_RNA"] = member["external_state_restricted_source_RNA"].astype(str).str.lower().eq("true")
external_members = (member.loc[member.external_state_restricted_source_RNA & member.signature_group.isin(candidate.signature_group)]
                    .groupby("signature_group")["secreted_protein_member"].apply(lambda x: ";".join(sorted(set(x)))).to_dict())
candidate["external_source_RNA_members"] = candidate.signature_group.map(external_members)
assert candidate.external_source_RNA_members.notna().all()

scores_by_cohort = {
    "GSE32225": ROOT / "results/82_bulk_primary_association_analysis/gse32225_tumor_bulk_scores.tsv",
    "GSE107943": ROOT / "results/84_gse107943_direction_and_clinical_sensitivity/gse107943_30_tumor_scores_and_metadata.tsv",
    "GSE26566": ROOT / "results/86_gse26566_third_bulk_directional_replication/gse26566_104_tumor_scores_and_metadata.tsv",
}
proxy32225 = pd.read_csv(ROOT / "results/83_bulk_stromal_endothelial_adjustment/gse32225_fibroblast_endothelial_proxy_scores.tsv", sep="\t", index_col=0)
fresh_matrices = {}
if (R3 / "external_inputs/GPL6104.annot.gz").exists():
    fresh_matrices["GSE107943"] = read_gse107943()
    fresh_matrices["GSE26566"] = read_gse26566()

rows = []
for cohort, path in scores_by_cohort.items():
    data = pd.read_csv(path, sep="\t", index_col=0)
    if cohort == "GSE32225":
        data = data.join(proxy32225, how="left")
    for rec in candidate.itertuples(index=False):
        genes = str(rec.external_source_RNA_members).split(";")
        if cohort in fresh_matrices:
            matrix, _ = fresh_matrices[cohort]
            s, present = score_matrix(matrix, genes)
            # The freshly reconstructed matrices are keyed by audited GSM accession.
            # Align by accession rather than assigning by row position, because the
            # legacy score tables may have a different order or a filtered subset.
            s = s.reindex(data.index)
        else:
            s, present = score(data, genes)
        outcome = "candidate_corrected"
        data[outcome] = s
        primary = fit(data, outcome, "VM_anchor_score")
        minus = fit(data, outcome, "VM_anchor_without_MMP2_score")
        rows.append({
            "cohort": cohort, "signature_group": rec.signature_group,
            "external_source_RNA_members": rec.external_source_RNA_members,
            "genes_present": ";".join(present), "n_requested_genes": len(genes), "n_present_genes": len(present),
            **{f"four_proxy_{k}": v for k, v in primary.items()},
            **{f"without_MMP2_four_proxy_{k}": v for k, v in minus.items()},
        })

results = pd.DataFrame(rows)
for cohort in results.cohort.unique():
    mask = results.cohort.eq(cohort)
    for p_col, fdr_col in [("four_proxy_p_value", "four_proxy_BH_FDR"), ("without_MMP2_four_proxy_p_value", "without_MMP2_four_proxy_BH_FDR")]:
        temp = results.loc[mask, [p_col]].copy()
        ok = temp[p_col].notna()
        vals = np.full(mask.sum(), np.nan)
        vals[ok.to_numpy()] = multipletests(temp.loc[ok, p_col], method="fdr_bh")[1]
        results.loc[mask, fdr_col] = vals

wide = results.pivot(index=["signature_group", "external_source_RNA_members"], columns="cohort", values="four_proxy_beta").reset_index()
wide["evaluable_all_three"] = wide[["GSE32225", "GSE107943", "GSE26566"]].notna().all(axis=1)
wide["positive_all_three"] = wide.evaluable_all_three & wide[["GSE32225", "GSE107943", "GSE26566"]].gt(0).all(axis=1)
wide["same_direction_all_three"] = wide.evaluable_all_three & ((wide[["GSE32225", "GSE107943", "GSE26566"]].gt(0).all(axis=1)) | (wide[["GSE32225", "GSE107943", "GSE26566"]].lt(0).all(axis=1)))

results.to_csv(OUT / "Table_U14_corrected_26_group_bulk_composition_audit.tsv", sep="\t", index=False)
wide.to_csv(OUT / "Table_U15_corrected_26_group_three_bulk_direction.tsv", sep="\t", index=False)
(OUT / "analysis_manifest.json").write_text(json.dumps({
    "source_of_truth": "Reviewer 3 corrected q80 external-compatible 26-group set",
    "testing_family": "All 26 corrected groups in each bulk cohort",
    "model": "candidate source-RNA-member score ~ continuous six-anchor score + pan-myeloid + fibroblast + endothelial expression proxies; HC3",
    "boundary": "Bulk covariance only; no cellular-source, secretion, communication or causal claim.",
}, indent=2), encoding="utf-8")
print(wide.to_string(index=False))
print("\nPositive in all three:", int(wide.positive_all_three.sum()))
