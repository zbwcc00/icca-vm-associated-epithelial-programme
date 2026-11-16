from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:\第六篇大论文")
REVISION = ROOT / "revision_round1_reviewer3_20261116"
INPUT = REVISION / "analysis_03_corrected_edger_secact"
OUT = REVISION / "analysis_04_corrected_secact_integration"
OUT.mkdir(parents=True, exist_ok=True)
MEMBERSHIP = ROOT / "results/62_secact_receiver_activity/secact_signature_group_membership.tsv"
OLD_CANDIDATES = ROOT / "results/63_secact_cross_cohort_receiver_candidates/secact_cross_cohort_positive_nonredundant_signature_groups.tsv"


def load_one(label: str) -> pd.DataFrame:
    zscore = pd.read_csv(INPUT / f"{label}_SecAct_zscore.tsv", sep="\t")
    pvalue = pd.read_csv(INPUT / f"{label}_SecAct_pvalue.tsv", sep="\t")
    merged = zscore.merge(pvalue, on="secreted_protein", validate="one_to_one")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    merged = merged.merge(membership, on="secreted_protein", validate="one_to_one")
    groups = (
        merged.groupby(["signature_group", "group_members", "group_size"], observed=True)
        .agg(zscore=("zscore", "first"), pvalue=("pvalue", "first"))
        .reset_index()
    )
    groups["FDR_657_groups"] = multipletests(groups["pvalue"], method="fdr_bh")[1]
    groups["analysis"] = label
    return groups


def main() -> None:
    labels = ["GSE138709_q80", "GSE138709_q90", "GSE125449_q80", "GSE125449_q90"]
    tables = {label: load_one(label) for label in labels}
    pd.concat(tables.values(), ignore_index=True).to_csv(
        OUT / "Table_R3_6_corrected_threshold_SecAct_657_group_statistics.tsv", sep="\t", index=False
    )
    threshold_correlations = []
    for cohort in ["GSE138709", "GSE125449"]:
        q80 = tables[f"{cohort}_q80"].set_index("signature_group")
        q90 = tables[f"{cohort}_q90"].set_index("signature_group")
        rho, p_value = stats.spearmanr(q80["zscore"], q90.loc[q80.index, "zscore"])
        threshold_correlations.append(
            {
                "cohort": cohort,
                "n_nonredundant_groups": len(q80),
                "spearman_zscore_q80_vs_q90": rho,
                "p_value": p_value,
                "positive_direction_both_thresholds": int(
                    ((q80["zscore"] > 0) & (q90.loc[q80.index, "zscore"] > 0)).sum()
                ),
                "FDR_lt_0_05_both_thresholds": int(
                    ((q80["FDR_657_groups"] < 0.05) & (q90.loc[q80.index, "FDR_657_groups"] < 0.05)).sum()
                ),
            }
        )
    correlation_table = pd.DataFrame(threshold_correlations)
    correlation_table.to_csv(
        OUT / "Table_R3_7_SecAct_q80_q90_threshold_stability.tsv", sep="\t", index=False
    )

    cross_rows = []
    candidate_sets: dict[int, set[str]] = {}
    for threshold in [80, 90]:
        discovery = tables[f"GSE138709_q{threshold}"].set_index("signature_group")
        external = tables[f"GSE125449_q{threshold}"].set_index("signature_group")
        common = discovery.index.intersection(external.index)
        merged = pd.DataFrame(
            {
                "signature_group": common,
                "discovery_zscore": discovery.loc[common, "zscore"].to_numpy(),
                "discovery_pvalue": discovery.loc[common, "pvalue"].to_numpy(),
                "discovery_FDR": discovery.loc[common, "FDR_657_groups"].to_numpy(),
                "external_zscore": external.loc[common, "zscore"].to_numpy(),
                "external_pvalue": external.loc[common, "pvalue"].to_numpy(),
                "external_FDR": external.loc[common, "FDR_657_groups"].to_numpy(),
                "group_members": discovery.loc[common, "group_members"].to_numpy(),
                "group_size": discovery.loc[common, "group_size"].to_numpy(),
            }
        )
        merged["positive_both"] = (merged["discovery_zscore"] > 0) & (merged["external_zscore"] > 0)
        merged["nominal_positive_both"] = merged["positive_both"] & (merged["discovery_pvalue"] < 0.05) & (merged["external_pvalue"] < 0.05)
        merged["global_FDR_positive_both"] = merged["positive_both"] & (merged["discovery_FDR"] < 0.05) & (merged["external_FDR"] < 0.05)
        merged["threshold"] = threshold / 100
        merged.to_csv(OUT / f"corrected_cross_cohort_SecAct_q{threshold}.tsv", sep="\t", index=False)
        candidate_sets[threshold] = set(merged.loc[merged["nominal_positive_both"], "signature_group"])
        cross_rows.append(
            {
                "threshold": threshold / 100,
                "nonredundant_groups": len(merged),
                "positive_both": int(merged["positive_both"].sum()),
                "nominal_positive_both": int(merged["nominal_positive_both"].sum()),
                "global_FDR_positive_both": int(merged["global_FDR_positive_both"].sum()),
            }
        )
    old = pd.read_csv(OLD_CANDIDATES, sep="\t")
    old_set = set(old["signature_group"].astype(str))
    overlap = candidate_sets[80] & candidate_sets[90]
    union = candidate_sets[80] | candidate_sets[90]
    comparison = pd.DataFrame(
        [
            {"comparison": "q80_vs_q90_corrected", "set_a": len(candidate_sets[80]), "set_b": len(candidate_sets[90]), "intersection": len(overlap), "union": len(union), "jaccard": len(overlap) / len(union) if union else 1.0},
            {"comparison": "old_vs_corrected_q80", "set_a": len(old_set), "set_b": len(candidate_sets[80]), "intersection": len(old_set & candidate_sets[80]), "union": len(old_set | candidate_sets[80]), "jaccard": len(old_set & candidate_sets[80]) / len(old_set | candidate_sets[80]) if old_set | candidate_sets[80] else 1.0},
            {"comparison": "old_vs_corrected_q90", "set_a": len(old_set), "set_b": len(candidate_sets[90]), "intersection": len(old_set & candidate_sets[90]), "union": len(old_set | candidate_sets[90]), "jaccard": len(old_set & candidate_sets[90]) / len(old_set | candidate_sets[90]) if old_set | candidate_sets[90] else 1.0},
        ]
    )
    pd.DataFrame(cross_rows).to_csv(
        OUT / "Table_R3_8_corrected_cross_cohort_SecAct_summary.tsv", sep="\t", index=False
    )
    comparison.to_csv(
        OUT / "Table_R3_9_corrected_SecAct_candidate_overlap.tsv", sep="\t", index=False
    )
    stable = pd.DataFrame({"signature_group": sorted(overlap)})
    stable["present_in_old_nominal_candidate_set"] = stable["signature_group"].isin(old_set)
    stable.to_csv(OUT / "corrected_q80_q90_stable_nominal_candidates.tsv", sep="\t", index=False)
    manifest = {
        "family": "657 nonredundant SecAct signature groups per cohort and threshold",
        "multiplicity": "Benjamini-Hochberg within each 657-group family",
        "operational_thresholds": [0.80, 0.90],
        "candidate_rule": "Positive z-score and nominal P<0.05 in both cohorts at the same operational threshold",
        "stable_candidate_rule": "Candidate rule satisfied at both 80th- and 90th-percentile definitions",
        "interpretation": "Internal threshold-sensitivity calibration only; not independent validation or causal communication evidence.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("Threshold correlations")
    print(correlation_table.to_string(index=False))
    print("\nCross-cohort summaries")
    print(pd.DataFrame(cross_rows).to_string(index=False))
    print("\nCandidate overlaps")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
