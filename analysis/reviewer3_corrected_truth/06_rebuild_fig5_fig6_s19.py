from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats


ROOT = Path(r"D:\第六篇大论文")
REVISION = ROOT / "revision_round1_reviewer3_20261116"
OUT = REVISION / "figures_revised"
SOURCE = REVISION / "figure_source_data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)
MULTI = ROOT / "results/104_oep001105_multioics_state_validation"
GENOMIC = ROOT / "results/105_oep001105_genomic_state_validation"
REPORTING = REVISION / "analysis_05_phosphopeptide_genomic_reporting"
SCORING = REVISION / "analysis_01_corrected_scoring_and_anchor_ablation"
SECACT = REVISION / "analysis_04_corrected_secact_integration"

COLORS = {
    "dark": "#263238",
    "grey": "#9AA5AE",
    "light_grey": "#E7ECEF",
    "vm": "#CC79A7",
    "blue": "#0072B2",
    "light_blue": "#B9DDF2",
    "teal": "#009E73",
    "orange": "#E69F00",
    "protein": "#D55E00",
    "red": "#A61B1B",
    "grid": "#E9ECEF",
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.titleweight": "bold",
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 300,
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
        }
    )


def save(figure: plt.Figure, name: str) -> None:
    figure.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(figure)


def panel_label(axis: plt.Axes, value: str, x: float = -0.12, y: float = 1.12) -> None:
    axis.text(x, y, value, transform=axis.transAxes, fontweight="bold", fontsize=10, va="top", clip_on=False)


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}".replace("e-0", "e-")
    return f"{value:.3f}"


def fig5() -> None:
    shared = pd.read_csv(MULTI / "shared_sample_vm_and_endothelial_scores.tsv", sep="\t").dropna(
        subset=["VM_RNA_log2TPM", "VM_protein_NA"]
    )
    anchors = pd.read_csv(MULTI / "vm_anchor_gene_cross_omic_correlations.tsv", sep="\t")
    clinical = pd.read_csv(MULTI / "oep001105_vm_clinical_association_tests.tsv", sep="\t")
    phospho = pd.read_csv(REPORTING / "Table_R3_10_RNA_phosphopeptide_correlations_with_CI.tsv", sep="\t")
    protein = anchors.query("left_layer == 'RNA_log2TPM' and right_layer == 'protein_NA'").copy()
    protein = protein.sort_values("rho")

    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.4), gridspec_kw={"hspace": 0.76, "wspace": 0.48})
    rho, p_value = stats.spearmanr(shared["VM_RNA_log2TPM"], shared["VM_protein_NA"])
    axes[0, 0].scatter(
        shared["VM_RNA_log2TPM"], shared["VM_protein_NA"], s=12, color=COLORS["vm"], alpha=0.62, linewidths=0
    )
    coefficients = np.polyfit(shared["VM_RNA_log2TPM"], shared["VM_protein_NA"], 1)
    values = np.linspace(shared["VM_RNA_log2TPM"].min(), shared["VM_RNA_log2TPM"].max(), 100)
    axes[0, 0].plot(values, np.polyval(coefficients, values), color=COLORS["dark"], linewidth=1.2)
    axes[0, 0].set(
        xlabel="Six-gene RNA programme score",
        ylabel="Matched six-protein score",
        title="Programme-level RNA–protein association",
    )
    axes[0, 0].text(
        0.01,
        -0.25,
        f"n={len(shared)}; Spearman ρ={rho:.3f}; P={format_p(p_value)}",
        transform=axes[0, 0].transAxes,
        fontsize=6.8,
        ha="left",
        va="top",
        clip_on=False,
    )

    axes[0, 1].barh(protein["gene"], protein["rho"], color=COLORS["protein"])
    axes[0, 1].set(xlabel="RNA–total-protein Spearman ρ", xlim=(0, 1), title="Anchor-level cross-layer correlations")
    axes[0, 1].grid(axis="x", color=COLORS["grid"])
    for position, row in enumerate(protein.itertuples(index=False)):
        axes[0, 1].text(row.rho + 0.02, position, f"n={int(row.n)}", va="center", fontsize=6.3)

    measured = phospho.query("availability == 'measured'").copy()
    positions = np.arange(len(measured))
    errors = np.vstack(
        [measured["rho"] - measured["rho_ci95_low"], measured["rho_ci95_high"] - measured["rho"]]
    )
    axes[1, 0].bar(positions, measured["rho"], color=COLORS["orange"], width=0.62)
    axes[1, 0].errorbar(positions, measured["rho"], yerr=errors, fmt="none", color=COLORS["dark"], capsize=3, linewidth=0.9)
    axes[1, 0].set(
        xticks=positions,
        xticklabels=measured["gene"],
        ylabel="RNA–phosphopeptide Spearman ρ",
        ylim=(0, 0.9),
        title="Measured phosphopeptide associations",
    )
    for position, row in enumerate(measured.itertuples(index=False)):
        axes[1, 0].text(
            position,
            0.03,
            f"n={int(row.n)}\nP={format_p(row.p_value)}\nBH={format_p(row.fdr_bh_within_pair)}",
            ha="center",
            va="bottom",
            fontsize=5.8,
            color=COLORS["dark"],
        )
    axes[1, 0].text(
        0.01,
        -0.25,
        "95% CIs use a Fisher-z approximation; CDH5, MMP2 and MMP9 were not measured.",
        transform=axes[1, 0].transAxes,
        fontsize=6.2,
        va="top",
    )

    clinical = clinical.loc[clinical["outcome"].eq("VM_RNA_log2TPM_255")].copy()
    clinical["display"] = (
        clinical["test"].str.replace("Spearman: ", "", regex=False).str.replace("Mann-Whitney: ", "", regex=False)
    )
    clinical = clinical.sort_values("fdr_bh_all_prespecified_tests")
    axes[1, 1].scatter(
        -np.log10(clinical["fdr_bh_all_prespecified_tests"]), np.arange(len(clinical)), s=30, color=COLORS["grey"]
    )
    axes[1, 1].axvline(-math.log10(0.05), color=COLORS["red"], linestyle="--", linewidth=0.9)
    axes[1, 1].set(
        yticks=np.arange(len(clinical)),
        yticklabels=clinical["display"],
        xlabel="−log10 BH-FDR",
        title="Available surgery-period covariates",
    )
    for axis, name in zip(axes.flat, "ABCD"):
        panel_label(axis, name)
    figure.suptitle(
        "Fig. 5 | Matched cross-layer associations of the six-gene VM-related programme",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.005,
        "Cross-layer correlation does not establish measurement agreement, phosphorylation occupancy or pathway activation.",
        ha="center",
        fontsize=6.8,
    )
    shared.to_csv(SOURCE / "Fig5_programme_scores.tsv", sep="\t", index=False)
    protein.to_csv(SOURCE / "Fig5_anchor_RNA_total_protein_correlations.tsv", sep="\t", index=False)
    phospho.to_csv(SOURCE / "Fig5_phosphopeptide_statistics.tsv", sep="\t", index=False)
    clinical.to_csv(SOURCE / "Fig5_clinical_tests.tsv", sep="\t", index=False)
    save(figure, "Fig5_cross_layer_associations_revised")


def fig6() -> None:
    raw = pd.read_csv(ROOT / "repository_release/source_data/Fig6_matched_genomic_context.tsv", sep="\t")
    displayed = pd.read_csv(REPORTING / "Fig6B_displayed_genomic_comparisons.tsv", sep="\t")
    burdens = pd.read_csv(REPORTING / "Table_R3_12_genomic_burden_correlations_with_CI.tsv", sep="\t")
    protein = pd.read_csv(GENOMIC / "prespecified_genomic_features_protein_directional_check.tsv", sep="\t")
    cnv_audit = pd.read_csv(GENOMIC / "genomewide_discrete_CNV_VM_RNA_association.tsv", sep="\t")

    figure = plt.figure(figsize=(11.2, 7.8))
    grid = figure.add_gridspec(2, 2, hspace=0.78, wspace=0.50)
    burden_axis = figure.add_subplot(grid[0, 0])
    event_axis = figure.add_subplot(grid[0, 1])
    protein_axis = figure.add_subplot(grid[1, 0])
    genome_axis = figure.add_subplot(grid[1, 1])

    burden_specs = [
        ("nonsynonymous_mutation_burden", COLORS["blue"], "Mutation burden", burdens.iloc[0]),
        ("non_diploid_CNV_fraction", COLORS["teal"], "Non-diploid CNV fraction", burdens.iloc[1]),
    ]
    for column, color, name, row in burden_specs:
        standardized = (raw[column] - raw[column].mean()) / raw[column].std(ddof=1)
        burden_axis.scatter(standardized, raw["VM_RNA_log2TPM"], s=9, color=color, alpha=0.48, label=name, linewidths=0)
        burden_axis.text(
            0.01,
            0.98 if name.startswith("Mutation") else 0.88,
            f"{name}: ρ={row.rho:.3f} [{row.rho_ci95_low:.3f}, {row.rho_ci95_high:.3f}], P={format_p(row.p_value)}",
            transform=burden_axis.transAxes,
            fontsize=5.9,
            va="top",
            color=color,
        )
    burden_axis.legend(frameon=False, fontsize=6.3, loc="lower left")
    burden_axis.set(
        xlabel="Standardized genomic burden",
        ylabel="Six-gene RNA programme score",
        title="Burden associations in matched RNA–WES/CNV tumours",
    )

    event_specs = [
        ("KRAS_mutant", "KRAS\nmutation", "nonsynonymous mutation: KRAS"),
        ("BAP1_mutant", "BAP1\nmutation", "nonsynonymous mutation: BAP1"),
        ("MYC_gain", "MYC\ngain", "MYC gain (GISTIC > 0)"),
        ("MCL1_gain", "MCL1\ngain", "MCL1 gain (GISTIC > 0)"),
    ]
    centers = []
    for index, (column, title, feature) in enumerate(event_specs):
        center = index * 2.3 + 1.0
        centers.append(center + 0.36)
        event = raw[column].astype(bool)
        reference = raw.loc[~event, "VM_RNA_log2TPM"]
        altered = raw.loc[event, "VM_RNA_log2TPM"]
        boxes = event_axis.boxplot(
            [reference, altered],
            positions=[center, center + 0.72],
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": COLORS["dark"], "linewidth": 1.1},
        )
        boxes["boxes"][0].set_facecolor(COLORS["light_grey"])
        boxes["boxes"][1].set_facecolor(COLORS["light_blue"])
        report = displayed.loc[displayed["feature"].eq(feature)].iloc[0]
        event_axis.text(center, 1.42, f"Ref\nn={len(reference)}", ha="center", va="bottom", fontsize=5.7)
        event_axis.text(center + 0.72, 1.42, f"Alt\nn={len(altered)}", ha="center", va="bottom", fontsize=5.7)
        event_axis.text(
            center + 0.36,
            1.84,
            f"Δ={report.effect_altered_minus_reference_median:.2f}\nCI [{report.ci95_low:.2f}, {report.ci95_high:.2f}]\nBH={format_p(report.family_BH_FDR)}",
            ha="center",
            va="bottom",
            fontsize=5.4,
        )
    event_axis.set(
        xticks=centers,
        xticklabels=[item[1] for item in event_specs],
        ylabel="Six-gene RNA programme score",
        ylim=(-1.45, 2.40),
        title="Predefined mutation/CNV associations",
    )
    event_axis.legend(
        handles=[Patch(facecolor=COLORS["light_grey"], label="Reference"), Patch(facecolor=COLORS["light_blue"], label="Altered")],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=2,
    )

    protein = protein.loc[protein["eligible"].astype(bool)].copy()
    protein["display"] = (
        protein["feature"].str.replace("nonsynonymous mutation: ", "", regex=False).replace(
            {
                "coding nonsynonymous mutation burden": "Mutation burden",
                "fraction of genes with non-diploid discrete CNV": "CNV burden",
                "MYC gain (GISTIC > 0)": "MYC gain",
                "MCL1 gain (GISTIC > 0)": "MCL1 gain",
            }
        )
    )
    positions = np.arange(len(protein))
    protein_axis.hlines(positions, 0, protein["median_difference_altered_minus_wildtype"], color="#CDD4D9")
    protein_axis.scatter(protein["median_difference_altered_minus_wildtype"], positions, s=38, color=COLORS["protein"])
    protein_axis.axvline(0, color=COLORS["dark"], linewidth=0.7)
    protein_axis.set(
        yticks=positions,
        yticklabels=protein["display"],
        xlabel="Matched protein score, altered − reference median",
        title="Protein-level directional checks",
    )
    protein_axis.text(
        0.01,
        -0.23,
        "Same-project cross-layer associations; not an independent cohort.",
        transform=protein_axis.transAxes,
        fontsize=6.2,
        va="top",
    )

    cnv_audit["minus_log10_p"] = -np.log10(cnv_audit["p_value"])
    cnv_audit["minus_log10_fdr"] = -np.log10(cnv_audit["fdr_bh_genomewide"])
    genome_axis.scatter(
        cnv_audit["minus_log10_p"], cnv_audit["minus_log10_fdr"], s=4, color=COLORS["grey"], alpha=0.45, linewidths=0
    )
    genome_axis.axhline(-math.log10(0.05), color=COLORS["red"], linestyle="--", linewidth=0.9)
    genome_axis.set(
        xlabel="−log10 nominal P",
        ylabel="−log10 genome-wide BH-FDR",
        title="Genome-wide discrete CNV audit, 22,845 genes",
    )
    genome_axis.text(0.02, 0.90, "0 genes at BH-FDR < 0.05", transform=genome_axis.transAxes, va="top", fontsize=6.4)

    for axis, name in zip([burden_axis, event_axis, protein_axis, genome_axis], "ABCD"):
        panel_label(axis, name)
    figure.suptitle(
        "Fig. 6 | Quantitative genomic context of the six-gene VM-related programme",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.005,
        "All genomic results are sample-level associations; they do not identify a genomic driver or causal VM mechanism.",
        ha="center",
        fontsize=6.8,
    )
    raw.to_csv(SOURCE / "Fig6_matched_genomic_context.tsv", sep="\t", index=False)
    displayed.to_csv(SOURCE / "Fig6_displayed_event_statistics.tsv", sep="\t", index=False)
    burdens.to_csv(SOURCE / "Fig6_burden_statistics.tsv", sep="\t", index=False)
    save(figure, "Fig6_quantitative_genomic_context_revised")


def fig_s19() -> None:
    threshold = pd.read_csv(
        SCORING / "Table_R3_1_corrected_threshold_and_leave_one_anchor_reclassification.tsv", sep="\t"
    )
    previous = pd.read_csv(SCORING / "Table_R3_2_previous_vs_corrected_full6_classification.tsv", sep="\t")
    held = pd.read_csv(SCORING / "Table_R3_3_held_out_anchor_patient_pseudobulk_tests.tsv", sep="\t")
    stability = pd.read_csv(SECACT / "Table_R3_7_SecAct_q80_q90_threshold_stability.tsv", sep="\t")
    overlap = pd.read_csv(SECACT / "Table_R3_9_corrected_SecAct_candidate_overlap.tsv", sep="\t")

    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.4), gridspec_kw={"hspace": 0.70, "wspace": 0.42})
    prevalence_rows = []
    for cohort in ["GSE138709", "GSE125449"]:
        old = previous.query("cohort == @cohort and classification == 'previous_reported_label'").iloc[0]
        q80 = threshold.query("cohort == @cohort and classification == 'full6' and VM_score_quantile == 0.8").iloc[0]
        q90 = threshold.query("cohort == @cohort and classification == 'full6' and VM_score_quantile == 0.9").iloc[0]
        prevalence_rows.extend(
            [
                {"cohort": cohort, "definition": "Previous", "fraction": old.high_fraction},
                {"cohort": cohort, "definition": "Corrected q80", "fraction": q80.high_fraction},
                {"cohort": cohort, "definition": "Corrected q90", "fraction": q90.high_fraction},
            ]
        )
    prevalence = pd.DataFrame(prevalence_rows)
    definitions = ["Previous", "Corrected q80", "Corrected q90"]
    x = np.arange(2)
    width = 0.23
    for index, definition in enumerate(definitions):
        values = prevalence.query("definition == @definition").set_index("cohort").loc[["GSE138709", "GSE125449"], "fraction"]
        axes[0, 0].bar(x + (index - 1) * width, values * 100, width, label=definition)
    axes[0, 0].set(xticks=x, xticklabels=["GSE138709", "GSE125449"], ylabel="Operational high-score cells (%)", title="Corrected six-gene upper-tail definitions")
    axes[0, 0].legend(frameon=False)

    ablation = threshold.query("classification != 'full6'").copy()
    matrix = ablation.pivot_table(index=["cohort", "VM_score_quantile"], columns="omitted_anchor", values="jaccard_vs_corrected_full6")
    matrix = matrix[["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]]
    image = axes[0, 1].imshow(matrix.to_numpy(), vmin=0, vmax=1, cmap="Blues", aspect="auto")
    axes[0, 1].set(
        xticks=np.arange(matrix.shape[1]),
        xticklabels=matrix.columns,
        yticks=np.arange(matrix.shape[0]),
        yticklabels=[f"{cohort} q{int(quantile*100)}" for cohort, quantile in matrix.index],
        title="De novo leave-one-anchor reclassification",
    )
    axes[0, 1].tick_params(axis="x", rotation=35)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            axes[0, 1].text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=6.2, color="white" if value > 0.65 else COLORS["dark"])
    figure.colorbar(image, ax=axes[0, 1], label="Jaccard versus full six-gene call", shrink=0.82)

    anchors = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
    y = np.arange(len(anchors))
    for offset, cohort, color in [(-0.12, "GSE138709", COLORS["blue"]), (0.12, "GSE125449", COLORS["teal"])]:
        values = held.query("cohort == @cohort").set_index("omitted_anchor").loc[anchors]
        axes[1, 0].scatter(values["median_high_minus_low_log1p_CPM"], y + offset, color=color, s=32, label=cohort)
    axes[1, 0].axvline(0, color=COLORS["dark"], linewidth=0.8)
    axes[1, 0].set(yticks=y, yticklabels=anchors, xlabel="Held-out gene median high − low log1p CPM", title="Held-out-anchor patient-pseudobulk audit")
    axes[1, 0].legend(frameon=False)
    axes[1, 0].text(0.01, -0.22, "No held-out test survived within-cohort BH correction.", transform=axes[1, 0].transAxes, fontsize=6.3, va="top")

    bars = axes[1, 1].bar(stability["cohort"], stability["spearman_zscore_q80_vs_q90"], color=[COLORS["blue"], COLORS["teal"]], width=0.55)
    axes[1, 1].set(ylim=(0, 1.05), ylabel="Spearman ρ", title="SecAct z-score stability across q80 and q90")
    for bar, value in zip(bars, stability["spearman_zscore_q80_vs_q90"]):
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{value:.3f}", ha="center", fontsize=7)
    overlap_row = overlap.query("comparison == 'q80_vs_q90_corrected'").iloc[0]
    axes[1, 1].text(
        0.5,
        0.10,
        f"Nominal cross-cohort candidates: q80={int(overlap_row.set_a)}, q90={int(overlap_row.set_b)}\nIntersection={int(overlap_row.intersection)}; Jaccard={overlap_row.jaccard:.3f}",
        transform=axes[1, 1].transAxes,
        ha="center",
        fontsize=6.5,
    )

    for axis, name in zip(axes.flat, "ABCD"):
        panel_label(axis, name)
    figure.suptitle(
        "Fig. S19 | Corrected six-gene scoring, anchor ablation and threshold sensitivity",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.005,
        "Upper-tail labels are operational summaries of a continuous high-score programme; they are not discrete cell states or morphological VM calls.",
        ha="center",
        fontsize=6.8,
    )
    prevalence.to_csv(SOURCE / "FigS19_prevalence.tsv", sep="\t", index=False)
    matrix.to_csv(SOURCE / "FigS19_anchor_ablation_jaccard.tsv", sep="\t")
    held.to_csv(SOURCE / "FigS19_held_out_anchor_tests.tsv", sep="\t", index=False)
    stability.to_csv(SOURCE / "FigS19_SecAct_threshold_stability.tsv", sep="\t", index=False)
    save(figure, "FigS19_corrected_scoring_anchor_ablation_threshold_sensitivity")


def main() -> None:
    style()
    fig5()
    fig6()
    fig_s19()
    print(f"Wrote revised figures to {OUT}")


if __name__ == "__main__":
    main()
