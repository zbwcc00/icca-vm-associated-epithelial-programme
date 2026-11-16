from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(r"D:\第六篇大论文")
R3 = ROOT / "revision_round1_reviewer3_20261116"
SCORING = R3 / "analysis_01_corrected_scoring_and_anchor_ablation"
OUT = R3 / "figures_revised"
SOURCE = R3 / "figure_source_data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
    }
)

COHORTS = {
    "GSE138709": SCORING / "GSE138709_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
    "GSE181878": ROOT
    / "results/29_gse181878_vm_triple_scoring/malignant_epithelial_vm_triple_scores.tsv.gz",
    "GSE125449": SCORING / "GSE125449_corrected_full6_and_leave_one_anchor_scores.tsv.gz",
}
LABELS = ["Module score", "UCell", "AUCell"]
COLORS = {"preliminary": "#56B4E9", "final": "#CC79A7", "excluded": "#B7C3CE"}


def locate_columns(frame: pd.DataFrame) -> tuple[list[str], list[str]]:
    corrected_scores = ["module_full6", "ucell_full6", "aucell_full6"]
    corrected_controls = ["module_endothelial", "ucell_endothelial", "aucell_endothelial"]
    legacy_scores = ["AddModuleScore_VM_anchor", "UCell_VM_anchor", "AUCell_VM_anchor"]
    legacy_controls = [
        "AddModuleScore_endothelial_control",
        "UCell_endothelial_negative_control",
        "AUCell_endothelial_negative_control",
    ]
    if all(column in frame for column in corrected_scores + corrected_controls):
        return corrected_scores, corrected_controls
    if all(column in frame for column in legacy_scores + legacy_controls):
        return legacy_scores, legacy_controls
    raise KeyError(f"Unrecognised scoring columns: {frame.columns.tolist()}")


frames: dict[str, pd.DataFrame] = {}
correlation_rows = []
frequency_rows = []
for cohort, path in COHORTS.items():
    frame = pd.read_csv(path, sep="\t")
    score_columns, control_columns = locate_columns(frame)
    renamed = frame[score_columns].copy()
    renamed.columns = LABELS
    frames[cohort] = renamed
    correlation = renamed.corr(method="spearman")
    for first in LABELS:
        for second in LABELS:
            correlation_rows.append(
                {"cohort": cohort, "score_1": first, "score_2": second, "spearman_rho": correlation.loc[first, second]}
            )

    score_thresholds = frame[score_columns].quantile(0.80)
    control_thresholds = frame[control_columns].quantile(0.90)
    preliminary = frame[score_columns].gt(score_thresholds).all(axis=1)
    endothelial_flag = frame[control_columns].gt(control_thresholds).any(axis=1)
    final = preliminary & ~endothelial_flag
    frequency_rows.extend(
        [
            {"cohort": cohort, "category": "Preliminary high-score", "n_cells": int(preliminary.sum()), "fraction": preliminary.mean()},
            {"cohort": cohort, "category": "Endothelial-control excluded", "n_cells": int((preliminary & endothelial_flag).sum()), "fraction": (preliminary & endothelial_flag).mean()},
            {"cohort": cohort, "category": "Final operational high-score", "n_cells": int(final.sum()), "fraction": final.mean()},
        ]
    )

correlations = pd.DataFrame(correlation_rows)
frequencies = pd.DataFrame(frequency_rows)
correlations.to_csv(SOURCE / "FigS2_corrected_score_implementation_correlations.tsv", sep="\t", index=False)
frequencies.to_csv(SOURCE / "FigS3_corrected_high_score_frequencies.tsv", sep="\t", index=False)

figure, axes = plt.subplots(1, 3, figsize=(8.2, 2.85), constrained_layout=True)
for axis, cohort in zip(axes, COHORTS):
    matrix = correlations.query("cohort == @cohort").pivot(index="score_1", columns="score_2", values="spearman_rho")
    matrix = matrix.loc[LABELS, LABELS]
    sns.heatmap(
        matrix,
        ax=axis,
        vmin=0,
        vmax=1,
        cmap=sns.light_palette("#0072B2", as_cmap=True),
        annot=True,
        fmt=".2f",
        square=True,
        linewidths=1,
        linecolor="white",
        cbar=False,
        annot_kws={"fontsize": 8},
    )
    axis.set_title(cohort)
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=35)
    axis.tick_params(axis="y", rotation=0)
figure.suptitle("Agreement among three implementations of the corrected six-gene score", fontsize=11, fontweight="bold")
figure.text(0.5, -0.02, "Agreement assesses implementation robustness; it is not independent biological validation.", ha="center", fontsize=8)
figure.savefig(OUT / "FigS2_corrected_score_implementation_agreement.pdf")
figure.savefig(OUT / "FigS2_corrected_score_implementation_agreement.png")
plt.close(figure)

figure, axis = plt.subplots(figsize=(6.6, 3.7), constrained_layout=True)
x = np.arange(len(COHORTS))
width = 0.24
categories = [
    ("Preliminary high-score", COLORS["preliminary"]),
    ("Endothelial-control excluded", COLORS["excluded"]),
    ("Final operational high-score", COLORS["final"]),
]
for offset, (category, color) in enumerate(categories):
    values = (
        frequencies.query("category == @category").set_index("cohort").loc[list(COHORTS), "fraction"].to_numpy() * 100
    )
    bars = axis.bar(x + (offset - 1) * width, values, width=width, color=color, label=category)
    for bar, value in zip(bars, values):
        axis.text(bar.get_x() + bar.get_width() / 2, value + 0.25, f"{value:.2f}%", ha="center", va="bottom", fontsize=7)
axis.set_xticks(x)
axis.set_xticklabels(list(COHORTS))
axis.set_ylabel("Cells among operational malignant epithelium (%)")
axis.set_title("Corrected q80 operational upper-tail frequencies")
axis.set_ylim(0, max(18, axis.get_ylim()[1] * 1.18))
axis.grid(axis="y", alpha=0.18)
axis.legend(loc="upper left")
axis.text(
    0.01,
    -0.22,
    "GSE138709 used all six anchors after MMP2, MMP9 and LAMC2 were restored from raw per-library UMI matrices.",
    transform=axis.transAxes,
    fontsize=8,
)
figure.savefig(OUT / "FigS3_corrected_operational_high_score_frequencies.pdf")
figure.savefig(OUT / "FigS3_corrected_operational_high_score_frequencies.png")
plt.close(figure)

print("Rebuilt corrected Supplementary Figures S2 and S3")


def panel_label(axis, label):
    axis.text(-0.12, 1.05, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top")


def save(figure, stem):
    figure.savefig(OUT / f"{stem}.pdf")
    figure.savefig(OUT / f"{stem}.png")
    plt.close(figure)


constraint = pd.read_csv(
    R3 / "analysis_04_corrected_secact_integration/Table_R3_13_corrected_q80_source_and_external_constraint.tsv",
    sep="\t",
)
q80_candidates = constraint.loc[constraint["nominal_positive_both"].astype(bool)].copy()
source_restricted = q80_candidates.loc[q80_candidates["has_source_restricted_member"].astype(bool)].copy()
external_count = int(source_restricted["external_RNA_compatible"].astype(bool).sum())
not_external_count = len(source_restricted) - external_count
figure, axis = plt.subplots(figsize=(6.4, 3.8))
labels = ["External RNA-compatible", "No external passing member"]
values = [external_count, not_external_count]
bars = axis.barh(
    labels,
    values,
    color=["#CC79A7", "#B7C3CE"],
    height=0.55,
)
for bar, value in zip(bars, values):
    axis.text(
        value + 0.5,
        bar.get_y() + bar.get_height() / 2,
        str(value),
        va="center",
        fontweight="bold",
    )
axis.set_xlabel("q80 source-restricted SecAct groups")
axis.set_xlim(0, max(values) * 1.22)
axis.invert_yaxis()
axis.set_title("External source-RNA compatibility audit", pad=10)
axis.text(
    0.5,
    -0.22,
    "GSE201425 provides an external source-RNA compatibility audit, not independent validation.",
    transform=axis.transAxes,
    fontsize=8,
    va="top",
    ha="center",
)
figure.subplots_adjust(left=0.34, right=0.97, bottom=0.28, top=0.84)
save(figure, "FigS5_corrected_external_source_RNA_audit")


threshold = pd.read_csv(
    SCORING / "Table_R3_1_corrected_threshold_and_leave_one_anchor_reclassification.tsv",
    sep="\t",
)
secact_summary = pd.read_csv(
    R3 / "analysis_04_corrected_secact_integration/Table_R3_8_corrected_cross_cohort_SecAct_summary.tsv",
    sep="\t",
)
stability = pd.read_csv(
    R3 / "analysis_04_corrected_secact_integration/Table_R3_7_SecAct_q80_q90_threshold_stability.tsv",
    sep="\t",
)
overlap = pd.read_csv(
    R3 / "analysis_04_corrected_secact_integration/Table_R3_9_corrected_SecAct_candidate_overlap.tsv",
    sep="\t",
)
r2 = ROOT / "revision_round1_reviewer2_20261115"
random_null = pd.read_csv(
    r2 / "analysis_03_detection_matched_random_benchmark/random_effect_distribution_all_cohorts.tsv.gz",
    sep="\t",
)
random_observed = pd.read_csv(
    r2 / "analysis_03_detection_matched_random_benchmark/Table_R2_6_detection_and_expression_matched_random_benchmark.tsv",
    sep="\t",
)
programme = pd.read_csv(
    r2 / "analysis_02_program_overlap/Table_R2_5a_VM_program_correlations_and_partial_correlations.tsv",
    sep="\t",
)

figure, axes = plt.subplots(2, 2, figsize=(7.4, 6.3), constrained_layout=True)
axis = axes[0, 0]
x = np.arange(2)
width = 0.23
for offset, column, label, color in [
    (-width, "nonredundant_groups", "All groups", "#B7C3CE"),
    (0, "nominal_positive_both", "Nominal positive", "#56B4E9"),
    (width, "global_FDR_positive_both", "Global BH-FDR positive", "#0072B2"),
]:
    bars = axis.bar(x + offset, secact_summary[column], width=width, color=color, label=label)
    for bar, value in zip(bars, secact_summary[column]):
        axis.text(bar.get_x() + bar.get_width() / 2, value + 10, str(int(value)), ha="center", fontsize=7)
axis.set_xticks(x, ["q80", "q90"])
axis.set_ylabel("SecAct signature groups")
axis.set_ylim(0, 730)
axis.legend(fontsize=7)
axis.set_title("Corrected global testing-family audit")
panel_label(axis, "A")

axis = axes[0, 1]
ablation = threshold.query("classification != 'full6'").copy()
for cohort, color, marker in [("GSE138709", "#0072B2", "o"), ("GSE125449", "#E69F00", "s")]:
    data = ablation.query("cohort == @cohort and VM_score_quantile == 0.8")
    axis.plot(data["omitted_anchor"], data["jaccard_vs_corrected_full6"], color=color, marker=marker, label=cohort)
axis.set_ylim(0.25, 1.03)
axis.set_ylabel("Jaccard vs full six-gene q80 label")
axis.tick_params(axis="x", rotation=35)
axis.legend(loc="lower left")
axis.set_title("De novo leave-one-anchor reclassification")
panel_label(axis, "B")

axis = axes[1, 0]
bulk = programme.query("dataset in ['GSE32225','GSE107943','GSE26566']").copy()
names = {
    "Hallmark_Hypoxia": "Hypoxia",
    "Hallmark_EMT": "EMT",
    "Hallmark_Angiogenesis": "Angio",
    "Hallmark_TGF_beta": "TGF-beta",
    "Puram_partial_EMT": "Puram pEMT",
}
matrix = bulk.pivot(index="dataset", columns="program", values="spearman_rho").reindex(columns=list(names))
matrix.columns = [names[column] for column in matrix.columns]
sns.heatmap(matrix, ax=axis, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", cbar=False, linewidths=0.6, linecolor="white")
axis.tick_params(axis="x", rotation=35)
axis.set_xlabel("")
axis.set_ylabel("")
axis.set_title("Published-programme calibration")
panel_label(axis, "C")

axis = axes[1, 1]
for cohort, color in [("GSE138709", "#0072B2"), ("GSE125449", "#E69F00")]:
    values = random_null.loc[random_null.cohort == cohort, "paired_Hedges_g"]
    sns.kdeplot(values, ax=axis, color=color, linewidth=1.8, label=cohort)
    row = random_observed.loc[random_observed.cohort == cohort].iloc[0]
    axis.axvline(row.six_anchor_paired_Hedges_g, color=color, linestyle="--", linewidth=1.4)
    axis.text(0.98, 0.92 if cohort == "GSE138709" else 0.78, f"{cohort}: P={row.empirical_two_sided_p:.3f}", transform=axis.transAxes, ha="right", color=color, fontsize=7)
axis.set_xlabel("Paired patient-pseudobulk Hedges g")
axis.set_ylabel("Density")
axis.legend(loc="upper left")
axis.set_title("Internal matched-random calibration")
panel_label(axis, "D")
figure.suptitle("Supplementary Figure S15 | Corrected robustness calibration", fontsize=11, fontweight="bold")
save(figure, "FigS15_corrected_robustness_calibration")


figure, axes = plt.subplots(2, 3, figsize=(10.5, 6.6), constrained_layout=True)
axis = axes[0, 0]
full = threshold.query("classification == 'full6'").copy()
for cohort, color in [("GSE138709", "#0072B2"), ("GSE125449", "#E69F00")]:
    data = full.query("cohort == @cohort")
    axis.plot(data.VM_score_quantile * 100, data.high_fraction * 100, marker="o", color=color, label=cohort)
    for row in data.itertuples():
        axis.text(row.VM_score_quantile * 100, row.high_fraction * 100 + 0.5, f"{int(row.n_patients_with_high)}/{int(row.n_patients)}p", ha="center", color=color, fontsize=7)
axis.set_xlabel("Programme-score percentile")
axis.set_ylabel("Operational high-score cells (%)")
axis.set_xticks([80, 90])
axis.legend()
axis.set_title("Threshold and patient coverage")
panel_label(axis, "A")

axis = axes[0, 1]
bars = axis.bar(stability.cohort, stability.spearman_zscore_q80_vs_q90, color=["#0072B2", "#E69F00"], width=0.55)
for bar, value in zip(bars, stability.spearman_zscore_q80_vs_q90):
    axis.text(bar.get_x() + bar.get_width() / 2, value - 0.08, f"ρ={value:.3f}", ha="center", color="white", fontweight="bold")
axis.set_ylim(0, 1.05)
axis.set_ylabel("SecAct z-score Spearman ρ")
axis.set_title("q80 vs q90 SecAct stability")
panel_label(axis, "B")

axis = axes[0, 2]
funnel = [657, 268, 58, 26]
labels = ["All", "Nominal\nboth", "Source\nmember", "External\nRNA"]
bars = axis.bar(labels, funnel, color=["#B7C3CE", "#56B4E9", "#009E73", "#CC79A7"])
for bar, value in zip(bars, funnel):
    axis.text(bar.get_x() + bar.get_width() / 2, value + 10, str(value), ha="center", fontsize=7)
axis.set_ylim(0, 730)
axis.set_ylabel("Groups")
axis.set_title("Corrected q80 evidence funnel")
panel_label(axis, "C")

axis = axes[1, 0]
row = overlap.query("comparison == 'q80_vs_q90_corrected'").iloc[0]
axis.barh(["q80", "q90", "Intersection"], [row.set_a, row.set_b, row.intersection], color=["#0072B2", "#E69F00", "#009E73"])
axis.set_xlabel("Nominal cross-cohort candidate groups")
axis.text(0.98, 0.10, f"Jaccard={row.jaccard:.3f}", transform=axis.transAxes, ha="right")
axis.set_title("Candidate-set threshold sensitivity")
panel_label(axis, "D")

axis = axes[1, 1]
axis.bar(["Source-restricted", "External compatible"], [58, 26], color=["#009E73", "#CC79A7"])
axis.set_ylabel("q80 groups")
axis.set_title("Source-RNA constraint")
for index, value in enumerate([58, 26]):
    axis.text(index, value + 1, str(value), ha="center")
panel_label(axis, "E")

axis = axes[1, 2]
compartment = pd.read_csv(
    ROOT / "revision_round1_reviewer1_20261115/analysis_outputs/PLAU_SERPINE1_compartment_expression_summary.tsv",
    sep="\t",
)
compartments = ["Myeloid", "Stromal", "Endothelial", "Malignant_epithelial_operational"]
matrix = compartment.query("compartment in @compartments").pivot(index="gene", columns="compartment", values="median_CPM").reindex(index=["PLAU", "SERPINE1"], columns=compartments)
sns.heatmap(np.log1p(matrix), ax=axis, cmap="YlGnBu", annot=matrix, fmt=".1f", cbar=False, linewidths=0.8, linecolor="white")
axis.set_xlabel("")
axis.set_ylabel("")
axis.tick_params(axis="x", rotation=35)
axis.set_title("Non-exclusive source RNA")
panel_label(axis, "F")
figure.suptitle("Supplementary Figure S16 | Corrected threshold and source-expression audits", fontsize=11, fontweight="bold")
save(figure, "FigS16_corrected_threshold_source_audits")

print("Rebuilt corrected Supplementary Figures S5, S15 and S16")
