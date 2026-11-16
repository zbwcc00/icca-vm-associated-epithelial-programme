from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"D:\第六篇大论文")
REVISION = ROOT / "revision_round1_reviewer3_20261116"
OUT = REVISION / "analysis_05_phosphopeptide_genomic_reporting"
OUT.mkdir(parents=True, exist_ok=True)
MULTI = ROOT / "results/104_oep001105_multioics_state_validation"
GENOMIC = ROOT / "results/105_oep001105_genomic_state_validation"
RAW_GENOMIC = ROOT / "repository_release/source_data/Fig6_matched_genomic_context.tsv"
SEED = 20261116
N_BOOT = 20000


def fisher_rho_ci(rho: float, n: int, confidence: float = 0.95) -> tuple[float, float]:
    if not np.isfinite(rho) or n <= 3 or abs(rho) >= 1:
        return np.nan, np.nan
    critical = stats.norm.ppf(1 - (1 - confidence) / 2)
    transformed = np.arctanh(rho)
    margin = critical / np.sqrt(n - 3)
    return float(np.tanh(transformed - margin)), float(np.tanh(transformed + margin))


def bootstrap_median_difference(reference: np.ndarray, altered: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    estimates = np.empty(N_BOOT, dtype=float)
    for index in range(N_BOOT):
        ref = rng.choice(reference, size=len(reference), replace=True)
        alt = rng.choice(altered, size=len(altered), replace=True)
        estimates[index] = np.median(alt) - np.median(ref)
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def phosphopeptide_table() -> pd.DataFrame:
    correlations = pd.read_csv(MULTI / "vm_anchor_gene_cross_omic_correlations.tsv", sep="\t")
    table = correlations.query("left_layer == 'RNA_log2TPM' and right_layer == 'phospho_NA'").copy()
    table = pd.DataFrame({"gene": ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]}).merge(
        table[["gene", "n", "rho", "p_value", "fdr_bh_within_pair"]], on="gene", how="left"
    )
    intervals = [fisher_rho_ci(row.rho, int(row.n)) if pd.notna(row.n) and row.n > 0 else (np.nan, np.nan) for row in table.itertuples()]
    table["rho_ci95_low"] = [item[0] for item in intervals]
    table["rho_ci95_high"] = [item[1] for item in intervals]
    table["availability"] = np.where(table["n"].fillna(0).gt(0), "measured", "not_available")
    table["ci_method"] = np.where(table["availability"].eq("measured"), "Fisher-z approximation", "not_applicable")
    return table


def complete_genomic_table() -> pd.DataFrame:
    mutation = pd.read_csv(GENOMIC / "prespecified_icca_driver_mutation_VM_RNA_associations.tsv", sep="\t")
    mutation = mutation.rename(
        columns={
            "median_difference_altered_minus_wildtype": "effect_altered_minus_reference_median",
            "fdr_bh_prespecified_genes": "family_BH_FDR",
        }
    )
    mutation["family"] = "prespecified_mutation_genes"
    mutation["n_reference"] = mutation["n"] - mutation["n_altered"]
    mutation["event_type"] = "mutation"
    mutation["ci95_low"] = np.nan
    mutation["ci95_high"] = np.nan

    cnv = pd.read_csv(GENOMIC / "prespecified_icca_driver_CNV_VM_RNA_associations.tsv", sep="\t")
    cnv = cnv.rename(
        columns={
            "median_difference_altered_minus_wildtype": "effect_altered_minus_reference_median",
            "fdr_bh_prespecified_events": "family_BH_FDR",
        }
    )
    cnv["family"] = "prespecified_CNV_events"
    cnv["n_reference"] = cnv["n"] - cnv["n_altered"]
    cnv["event_type"] = "copy_number"
    cnv["ci95_low"] = np.nan
    cnv["ci95_high"] = np.nan

    raw = pd.read_csv(RAW_GENOMIC, sep="\t")
    raw_map = {
        "nonsynonymous mutation: KRAS": "KRAS_mutant",
        "nonsynonymous mutation: BAP1": "BAP1_mutant",
        "MYC gain (GISTIC > 0)": "MYC_gain",
        "MCL1 gain (GISTIC > 0)": "MCL1_gain",
    }
    rng = np.random.default_rng(SEED)
    for frame in [mutation, cnv]:
        for row_index, feature in frame["feature"].items():
            if feature not in raw_map:
                continue
            event = raw[raw_map[feature]].astype(bool)
            reference = raw.loc[~event, "VM_RNA_log2TPM"].dropna().to_numpy(float)
            altered = raw.loc[event, "VM_RNA_log2TPM"].dropna().to_numpy(float)
            low, high = bootstrap_median_difference(reference, altered, rng)
            frame.loc[row_index, "ci95_low"] = low
            frame.loc[row_index, "ci95_high"] = high
    columns = [
        "family", "event_type", "feature", "gene", "event", "n", "n_reference", "n_altered",
        "effect_altered_minus_reference_median", "ci95_low", "ci95_high", "p_value", "family_BH_FDR",
        "eligible", "mann_whitney_u",
    ]
    for frame in [mutation, cnv]:
        for column in columns:
            if column not in frame:
                frame[column] = np.nan
    return pd.concat([mutation[columns], cnv[columns]], ignore_index=True)


def burden_table() -> pd.DataFrame:
    mutation = pd.read_csv(GENOMIC / "mutation_burden_VM_RNA_association.tsv", sep="\t")
    cnv = pd.read_csv(GENOMIC / "cnv_burden_VM_RNA_association.tsv", sep="\t")
    output = pd.concat([mutation, cnv], ignore_index=True)
    intervals = [fisher_rho_ci(float(row.rho), int(row.n)) for row in output.itertuples()]
    output["rho_ci95_low"] = [item[0] for item in intervals]
    output["rho_ci95_high"] = [item[1] for item in intervals]
    output["ci_method"] = "Fisher-z approximation"
    return output


def main() -> None:
    phospho = phosphopeptide_table()
    genomic = complete_genomic_table()
    burdens = burden_table()
    phospho.to_csv(OUT / "Table_R3_10_RNA_phosphopeptide_correlations_with_CI.tsv", sep="\t", index=False)
    genomic.to_csv(OUT / "Table_R3_11_complete_prespecified_genomic_comparisons.tsv", sep="\t", index=False)
    burdens.to_csv(OUT / "Table_R3_12_genomic_burden_correlations_with_CI.tsv", sep="\t", index=False)
    displayed = genomic[genomic["feature"].isin([
        "nonsynonymous mutation: KRAS",
        "nonsynonymous mutation: BAP1",
        "MYC gain (GISTIC > 0)",
        "MCL1 gain (GISTIC > 0)",
    ])].copy()
    displayed.to_csv(OUT / "Fig6B_displayed_genomic_comparisons.tsv", sep="\t", index=False)
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "phosphopeptide_CI": "Approximate 95% CI for Spearman rho using Fisher z transformation.",
                "displayed_genomic_event_CI": f"Percentile bootstrap CI for altered-minus-reference median difference; {N_BOOT} resamples; seed={SEED}.",
                "burden_CI": "Approximate 95% CI for Spearman rho using Fisher z transformation.",
                "multiplicity": "BH values are retained within the original prespecified mutation-gene family or CNV-event family.",
                "genomewide_CNV": "Separate 22,845-gene family reported elsewhere; no BH-FDR-positive genes.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print("Phosphopeptide reporting")
    print(phospho.to_string(index=False))
    print("\nDisplayed genomic comparisons")
    print(displayed.to_string(index=False))
    print("\nBurden reporting")
    print(burdens.to_string(index=False))


if __name__ == "__main__":
    main()
