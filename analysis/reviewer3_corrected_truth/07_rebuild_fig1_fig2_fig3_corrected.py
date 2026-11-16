from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch


ROOT = Path(r"D:\第六篇大论文")
R3 = ROOT / "revision_round1_reviewer3_20261116"
RESULTS = ROOT / "results"
REPOSITORY_SOURCE = ROOT / "repository_release" / "source_data"
OUT = R3 / "figures_revised"
SOURCE = R3 / "figure_source_data"
CORRECTED = R3 / "analysis_01_corrected_scoring_and_anchor_ablation"
SECACT = R3 / "analysis_04_corrected_secact_integration"

COLORS = {
    "dark": "#263238",
    "grey": "#6C757D",
    "light_grey": "#E9ECEF",
    "vm": "#CC79A7",
    "blue": "#0072B2",
    "teal": "#009E73",
    "orange": "#E69F00",
    "protein": "#D55E00",
    "red": "#A61B1B",
    "pale_pink": "#F9EEF5",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 300,
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(figure: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(figure)


def panel_label(axis: plt.Axes, label: str, x: float = -0.12, y: float = 1.12) -> None:
    axis.text(
        x,
        y,
        label,
        transform=axis.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.1},
        clip_on=False,
    )


def figure_1() -> None:
    figure, axis = plt.subplots(figsize=(7.7, 5.1))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.add_patch(
        FancyBboxPatch(
            (0.015, 0.015),
            0.97,
            0.97,
            boxstyle="round,pad=0.01,rounding_size=0.02",
            linewidth=0.6,
            facecolor="#FCFCFD",
            edgecolor="#D9DEE5",
        )
    )
    axis.text(
        0.04,
        0.945,
        "Fig. 1 | Evidence architecture of a six-gene VM-related epithelial programme",
        fontsize=8.8,
        fontweight="bold",
        va="center",
    )
    axis.text(
        0.04,
        0.912,
        "Public iCCA cohorts are integrated in a patient-aware framework; arrows show evidence flow, not causality.",
        fontsize=7,
        color=COLORS["grey"],
        va="center",
    )
    axis.plot([0.04, 0.96], [0.885, 0.885], color="#D9DEE5", linewidth=0.8)

    def card(x, y, width, height, title, body, accent, dashed=False):
        axis.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.018",
                linewidth=0.9,
                linestyle="--" if dashed else "-",
                facecolor="#FFFFFF",
                edgecolor=accent,
                zorder=2,
            )
        )
        axis.add_patch(
            FancyBboxPatch(
                (x, y + height - 0.043),
                width,
                0.043,
                boxstyle="round,pad=0.012,rounding_size=0.018",
                linewidth=0,
                facecolor=accent,
                edgecolor=accent,
                zorder=3,
            )
        )
        axis.add_patch(Circle((x + 0.022, y + height - 0.0215), 0.008, facecolor="#FFFFFF", edgecolor="none", zorder=4))
        axis.text(x + 0.04, y + height - 0.0215, title, fontsize=6.2, fontweight="bold", color="#FFFFFF", va="center", zorder=4)
        axis.text(x + width / 2, y + (height - 0.043) / 2, body, fontsize=6.1, ha="center", va="center", linespacing=1.35, zorder=4)

    card(0.045, 0.62, 0.245, 0.19, "01  SINGLE-CELL MAP", "3 iCCA scRNA-seq cohorts\n15 patients; malignant epithelium\nupper-tail programme coverage", COLORS["blue"])
    card(0.71, 0.62, 0.245, 0.19, "02  BULK CONSISTENCY", "3 independent bulk cohorts\n6 directionally consistent groups\ncomposition-aware sensitivity", COLORS["teal"])
    card(0.045, 0.305, 0.245, 0.17, "03  HYPOTHESIS CONTEXT", "Receiver activity plus\nexternal source-RNA compatibility\nmyeloid candidate context", COLORS["orange"], dashed=True)
    card(0.71, 0.305, 0.245, 0.17, "04  MATCHED OMICS", "OEP001105 iCCA tumours\nRNA–protein: n=208; ρ=0.709\nRNA–WES/CNV: n=249", COLORS["protein"])

    axis.add_patch(Ellipse((0.5, 0.57), 0.35, 0.25, facecolor="#F9D8E8", edgecolor=COLORS["vm"], linewidth=1.4, zorder=3))
    axis.text(0.5, 0.627, "VM-related malignant epithelial\nhigh-score programme", ha="center", va="center", fontsize=8.3, fontweight="bold", color="#5D254C", zorder=4)
    axis.text(0.5, 0.552, "six literature-derived anchors", ha="center", va="center", fontsize=6.5, color="#5D254C", zorder=4)
    for x, label in zip([0.425, 0.5, 0.575], ["CDH5", "EPHA2", "MCAM"]):
        axis.add_patch(FancyBboxPatch((x - 0.031, 0.50), 0.062, 0.027, boxstyle="round,pad=0.004,rounding_size=0.012", linewidth=0.3, facecolor="#FFFFFF", edgecolor="#D89ABB", zorder=4))
        axis.text(x, 0.5135, label, ha="center", va="center", fontsize=5.6, color="#5D254C", zorder=5)
    for x, label in zip([0.425, 0.5, 0.575], ["MMP2", "MMP9", "LAMC2"]):
        axis.add_patch(FancyBboxPatch((x - 0.031, 0.465), 0.062, 0.027, boxstyle="round,pad=0.004,rounding_size=0.012", linewidth=0.3, facecolor="#FFFFFF", edgecolor="#D89ABB", zorder=4))
        axis.text(x, 0.4785, label, ha="center", va="center", fontsize=5.6, color="#5D254C", zorder=5)

    def flow(start, end, color, linestyle="-"):
        axis.add_patch(FancyArrowPatch(start, end, connectionstyle="arc3,rad=0.12", arrowstyle="Simple,tail_width=0.35,head_width=4,head_length=5", linewidth=0.5, linestyle=linestyle, color=color, alpha=0.85, zorder=1))

    flow((0.29, 0.69), (0.355, 0.64), COLORS["blue"])
    flow((0.71, 0.69), (0.645, 0.64), COLORS["teal"])
    flow((0.29, 0.385), (0.36, 0.48), COLORS["orange"], "--")
    flow((0.71, 0.385), (0.64, 0.48), COLORS["protein"])
    axis.text(0.5, 0.395, "Cross-scale integration", ha="center", va="center", fontsize=7, color=COLORS["grey"])

    axis.add_patch(FancyBboxPatch((0.045, 0.095), 0.91, 0.12, boxstyle="round,pad=0.014,rounding_size=0.018", linewidth=0, facecolor="#EEF1F4", zorder=2))
    axis.text(0.065, 0.172, "FRAMEWORK-LEVEL CONCLUSION", fontsize=7, fontweight="bold", color=COLORS["dark"], va="center")
    axis.text(0.065, 0.145, "A continuous VM-related epithelial programme is measurable across heterogeneous iCCA cohorts.", fontsize=6.4, color=COLORS["dark"], va="center")
    axis.text(0.065, 0.118, "Within-project RNA–protein correlation; myeloid and genomic association-level context.", fontsize=6.4, color=COLORS["dark"], va="center")
    axis.text(0.5, 0.045, "Boundary: no vasculogenic-mimicry channels/perfusion, secretion/communication, or phosphorylation activation;\nno treatment resistance, independent prognosis, or drug vulnerability.", ha="center", va="center", fontsize=5.65, color="#9B1C1C")
    save_figure(figure, "Fig1_study_design_and_claim_boundaries_revised")


def figure_2() -> None:
    old = pd.read_csv(REPOSITORY_SOURCE / "Fig2_reconstructed_cell_scores_and_embeddings.tsv.gz", sep="\t")
    cohort_data = {}
    summaries = []
    patient_tables = []
    for cohort in ["GSE138709", "GSE181878", "GSE125449"]:
        base = old.loc[old["cohort"].eq(cohort)].copy()
        if cohort in {"GSE138709", "GSE125449"}:
            corrected = pd.read_csv(CORRECTED / f"{cohort}_corrected_full6_and_leave_one_anchor_scores.tsv.gz", sep="\t")
            keep = corrected[["cell_id", "patient_id", "module_full6", "ucell_full6", "aucell_full6", "corrected_full6_high_q80"]]
            base = base.drop(columns=["patient_id"], errors="ignore").merge(keep, on="cell_id", how="inner", validate="one_to_one")
            base["high_score_q80"] = base["corrected_full6_high_q80"].astype(bool)
            base["module_display"] = base["module_full6"]
            base["ucell_display"] = base["ucell_full6"]
            base["aucell_display"] = base["aucell_full6"]
        else:
            base["high_score_q80"] = base["VM_high_candidate"].astype(bool)
            base["module_display"] = base["AddModuleScore_VM_anchor"]
            base["ucell_display"] = base["UCell_VM_anchor"]
            base["aucell_display"] = base["AUCell_VM_anchor"]
        cohort_data[cohort] = base
        summaries.append(
            {
                "cohort": cohort,
                "patients": base["patient_id"].nunique(),
                "malignant_cells": len(base),
                "high_score_q80_cells": int(base["high_score_q80"].sum()),
                "high_score_q80_fraction": float(base["high_score_q80"].mean()),
                "patients_with_high_score": int(base.loc[base["high_score_q80"], "patient_id"].nunique()),
            }
        )

    summary = pd.DataFrame(summaries)
    figure = plt.figure(figsize=(11.0, 9.4))
    grid = figure.add_gridspec(4, 3, height_ratios=[0.58, 1.08, 1.06, 0.98], hspace=0.70, wspace=0.42)
    flow_axis = figure.add_subplot(grid[0, :])
    flow_axis.axis("off")
    flow_axis.set_xlim(0, 1)
    flow_axis.set_ylim(0, 1)
    flow_axis.text(0.0, 0.98, "A", fontsize=10, fontweight="bold", va="top")
    flow_axis.text(0.07, 0.98, "Cohort-aware operational upper-tail definition", fontsize=8.5, fontweight="bold", va="top")
    for position, row in zip([0.07, 0.38, 0.69], summary.itertuples(index=False)):
        flow_axis.add_patch(FancyBboxPatch((position, 0.18), 0.24, 0.56, boxstyle="round,pad=0.018,rounding_size=0.03", facecolor=COLORS["pale_pink"], edgecolor="#D7AEC6", linewidth=0.7))
        flow_axis.text(position + 0.12, 0.63, row.cohort, ha="center", fontsize=8, fontweight="bold")
        flow_axis.text(position + 0.12, 0.45, f"{row.patients} patient{'s' if row.patients != 1 else ''} | {row.malignant_cells:,} cells", ha="center", fontsize=6.6)
        flow_axis.text(position + 0.12, 0.31, f"{row.high_score_q80_cells:,} q80 high-score ({row.high_score_q80_fraction * 100:.1f}%)", ha="center", fontsize=6.4, color=COLORS["dark"])
    flow_axis.text(0.50, 0.01, "Operational malignant epithelium → three correlated score implementations → endothelial-control exclusion", ha="center", fontsize=6.9, color=COLORS["dark"])

    embedding_axes = [figure.add_subplot(grid[1, column]) for column in range(3)]
    score_axes = [figure.add_subplot(grid[2, column]) for column in range(3)]
    coverage_axes = [figure.add_subplot(grid[3, column]) for column in range(3)]
    for axis, (cohort, data) in zip(embedding_axes, cohort_data.items()):
        background = data.loc[~data["high_score_q80"]]
        high = data.loc[data["high_score_q80"]]
        axis.scatter(background["UMAP_1"], background["UMAP_2"], s=1.8, c="#C9D0D6", alpha=0.45, linewidths=0)
        axis.scatter(high["UMAP_1"], high["UMAP_2"], s=3.8, c=COLORS["vm"], alpha=0.78, linewidths=0)
        axis.set_title(cohort, pad=3)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xlabel("UMAP 1", labelpad=1)
        axis.set_ylabel("UMAP 2", labelpad=1)
    for axis, (cohort, data) in zip(score_axes, cohort_data.items()):
        background = data.loc[~data["high_score_q80"]]
        high = data.loc[data["high_score_q80"]]
        axis.scatter(background["ucell_display"], background["aucell_display"], s=2.4, c="#BFC8CF", alpha=0.35, linewidths=0)
        axis.scatter(high["ucell_display"], high["aucell_display"], s=4.2, c=COLORS["vm"], alpha=0.75, linewidths=0)
        axis.set_title(cohort, pad=3)
        axis.set_xlabel("UCell six-gene score")
        axis.set_ylabel("AUCell six-gene score")
        axis.grid(color=COLORS["light_grey"], linewidth=0.5)
    for axis, (cohort, data) in zip(coverage_axes, cohort_data.items()):
        coverage = data.groupby("patient_id", as_index=False).agg(malignant_cells=("cell_id", "size"), high_score_cells=("high_score_q80", "sum"))
        coverage["high_score_fraction"] = coverage["high_score_cells"] / coverage["malignant_cells"]
        coverage["cohort"] = cohort
        coverage = coverage.sort_values("high_score_fraction").reset_index(drop=True)
        patient_tables.append(coverage)
        positions = np.arange(len(coverage))
        sizes = 14 + 75 * np.sqrt(coverage["malignant_cells"] / coverage["malignant_cells"].max())
        axis.scatter(coverage["high_score_fraction"] * 100, positions, s=sizes, c=COLORS["vm"], edgecolor="white", linewidth=0.35, alpha=0.85)
        axis.set_yticks(positions)
        axis.set_yticklabels(coverage["patient_id"].astype(str).str.replace("GSE125449_", "", regex=False), fontsize=5.8)
        axis.set_title(cohort, pad=3)
        axis.set_xlabel("Operational q80 high-score cells (%)")
        axis.set_xlim(left=0)
        axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    panel_label(embedding_axes[0], "B")
    panel_label(score_axes[0], "C")
    panel_label(coverage_axes[0], "D")
    figure.suptitle("Fig. 2 | Cross-cohort distribution and patient coverage of a six-gene VM-related programme", x=0.01, y=0.995, ha="left", fontsize=11, fontweight="bold")
    figure.text(0.5, 0.012, "Purple cells meet the operational q80 rule after endothelial-control exclusion. The upper-tail label summarizes a continuous programme and is not an independent biological state or morphological VM call.", ha="center", fontsize=6.7)
    summary.to_csv(SOURCE / "Fig2_corrected_cohort_summary.tsv", sep="\t", index=False)
    pd.concat(patient_tables, ignore_index=True).to_csv(SOURCE / "Fig2_corrected_patient_coverage.tsv", sep="\t", index=False)
    pd.concat([data.assign(cohort=cohort) for cohort, data in cohort_data.items()], ignore_index=True).to_csv(SOURCE / "Fig2_corrected_cell_scores_and_embeddings.tsv.gz", sep="\t", index=False, compression="gzip")
    save_figure(figure, "Fig2_corrected_high_score_programme")


def h5ad_embedding(path: Path, columns: list[str]) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    coordinates = np.asarray(adata.obsm["X_umap"])
    data = adata.obs.loc[:, columns].copy()
    data["cell_id"] = data.index.astype(str)
    data["UMAP_1"] = coordinates[:, 0]
    data["UMAP_2"] = coordinates[:, 1]
    adata.file.close()
    return data.reset_index(drop=True)


def categorical_embedding(axis, data, column, title, palette):
    for category, color in palette.items():
        selected = data.loc[data[column].astype(str).eq(category)]
        if selected.empty:
            continue
        axis.scatter(selected["UMAP_1"], selected["UMAP_2"], s=3.0, color=color, alpha=0.72, linewidths=0, label=category)
    axis.set_title(title, pad=3)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlabel("UMAP 1", labelpad=1)
    axis.set_ylabel("UMAP 2", labelpad=1)


def figure_3() -> None:
    myeloid_data = h5ad_embedding(RESULTS / "61_myeloid_stroma_patient_coverage_audit" / "gse189903_myeloid_fine_state_audited.h5ad", ["fine_state", "case_prefix"])
    external_data = h5ad_embedding(RESULTS / "76_gse201425_independent_myeloid_identity_audit" / "gse201425_icca_myeloid_reference_mapped.h5ad", ["reference_state_call", "patient_id"])
    myeloid_coverage = pd.read_csv(RESULTS / "61_myeloid_stroma_patient_coverage_audit" / "myeloid_fine_state_patient_coverage.tsv", sep="\t")
    external_support = pd.read_csv(RESULTS / "76_gse201425_independent_myeloid_identity_audit" / "external_identity_support_summary.tsv", sep="\t")
    corrected = pd.read_csv(SECACT / "corrected_cross_cohort_SecAct_q80.tsv", sep="\t")
    source_groups = pd.read_csv(RESULTS / "66_secact_state_restricted_myeloid_prioritization" / "secact_groups_with_state_restricted_myeloid_members.tsv", sep="\t")
    source_members = pd.read_csv(RESULTS / "78_secact_multilayer_evidence_matrix" / "secact_source_member_multilayer_evidence_matrix.tsv", sep="\t")
    external_group = source_members.groupby("signature_group", as_index=False).agg(external_RNA_compatible=("external_state_restricted_source_RNA", "max"))
    source_lookup = source_groups[["signature_group", "state_restricted_source_members", "state_restricted_myeloid_states"]].drop_duplicates("signature_group")
    constrained = corrected.merge(source_lookup, on="signature_group", how="left").merge(external_group, on="signature_group", how="left")
    constrained["minimum_zscore"] = constrained[["discovery_zscore", "external_zscore"]].min(axis=1)
    constrained["mean_zscore"] = constrained[["discovery_zscore", "external_zscore"]].mean(axis=1)
    constrained["has_source_restricted_member"] = constrained["state_restricted_source_members"].notna()
    constrained["external_RNA_compatible"] = constrained["external_RNA_compatible"].fillna(False).astype(bool)
    candidates = constrained.loc[constrained["nominal_positive_both"].astype(bool)].copy()
    source_candidates = candidates.loc[candidates["has_source_restricted_member"]].copy()
    external_candidates = source_candidates.loc[source_candidates["external_RNA_compatible"]].copy()
    funnel_values = [657, len(candidates), len(source_candidates), len(external_candidates)]
    funnel_labels = ["Nonredundant SecAct groups", "Nominal positive in both scRNA cohorts", "State-restricted myeloid-source member", "External patient×state RNA-compatible"]

    eligible_states = myeloid_coverage.loc[myeloid_coverage["eligible_for_state_specific_secact"].astype(bool)].copy()
    myeloid_palette = {
        "C1QC_TAM_stress_associated": "#8E5EA2",
        "FCGR3A_macrophage_like": "#4C78A8",
        "HLAII_APC_like": "#72B7B2",
        "IL1B_CCL3_inflammatory_monocyte_like": "#E45756",
        "S100A8_S100A9_monocyte_like": "#F2CF5B",
        "HLAII_APC_ribosomal_like": "#B9B9B9",
        "T_cell_contamination_like": "#D6D6D6",
    }
    external_palette = {
        "C1QC_TAM_stress_associated": "#8E5EA2",
        "FCGR3A_macrophage_like": "#4C78A8",
        "HLAII_APC_like": "#72B7B2",
        "IL1B_CCL3_inflammatory_monocyte_like": "#E45756",
        "indeterminate_reference_transfer": "#D0D5D9",
    }
    figure = plt.figure(figsize=(11.0, 10.2))
    grid = figure.add_gridspec(3, 2, height_ratios=[1.0, 0.98, 1.10], hspace=0.62, wspace=0.40)
    myeloid_axis = figure.add_subplot(grid[0, 0])
    coverage_axis = figure.add_subplot(grid[0, 1])
    funnel_axis = figure.add_subplot(grid[1, 0])
    external_axis = figure.add_subplot(grid[1, 1])
    lollipop_axis = figure.add_subplot(grid[2, :])

    categorical_embedding(myeloid_axis, myeloid_data, "fine_state", "GSE189903 myeloid fine states", myeloid_palette)
    handles, labels = myeloid_axis.get_legend_handles_labels()
    eligible = set(eligible_states["fine_state"].astype(str))
    myeloid_axis.legend([h for h, l in zip(handles, labels) if l in eligible], [l.replace("_", " ") for l in labels if l in eligible], frameon=False, markerscale=2.0, fontsize=4.7, loc="upper center", bbox_to_anchor=(0.50, -0.19), ncol=2)
    states = eligible_states["fine_state"].astype(str).tolist()
    heatmap = eligible_states.loc[:, ["cells_1C", "cells_2C", "cells_3C"]].to_numpy(dtype=float)
    image = coverage_axis.imshow(heatmap, cmap="Purples", aspect="auto")
    coverage_axis.set_xticks(range(3))
    coverage_axis.set_xticklabels(["Patient 1C", "Patient 2C", "Patient 3C"])
    coverage_axis.set_yticks(range(len(states)))
    coverage_axis.set_yticklabels([state.replace("_", " ") for state in states], fontsize=5.8)
    coverage_axis.set_title("Patient coverage of eligible reference states", pad=3)
    for row in range(heatmap.shape[0]):
        for column in range(heatmap.shape[1]):
            coverage_axis.text(column, row, f"{int(heatmap[row, column])}", ha="center", va="center", fontsize=6.2, color="white" if heatmap[row, column] > np.nanmedian(heatmap) else COLORS["dark"])
    coverage_axis.figure.colorbar(image, ax=coverage_axis, fraction=0.045, pad=0.02, label="Cells")

    funnel_axis.barh(np.arange(4)[::-1], funnel_values, color=["#B7C3CE", COLORS["blue"], COLORS["teal"], COLORS["vm"]], height=0.62)
    funnel_axis.set_yticks(np.arange(4)[::-1])
    funnel_axis.set_yticklabels(funnel_labels, fontsize=6.1)
    funnel_axis.set_xlabel("Number of groups")
    funnel_axis.set_title("Corrected q80 evidence-constraining workflow", pad=3)
    for display_index, value in enumerate(funnel_values[::-1]):
        funnel_axis.text(value + max(funnel_values) * 0.012, display_index, str(value), va="center", fontsize=6.8, fontweight="bold")
    funnel_axis.set_xlim(0, max(funnel_values) * 1.14)
    funnel_axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)

    categorical_embedding(external_axis, external_data, "reference_state_call", "GSE201425 external reference transfer", external_palette)
    support_text = external_support.assign(state=external_support["state"].str.replace("_", " ", regex=False)).loc[:, ["state", "external_identity_support"]]
    external_axis.text(1.02, 0.98, "External identity support", transform=external_axis.transAxes, va="top", fontsize=6.2, fontweight="bold")
    for display_index, support_row in enumerate(support_text.itertuples(index=False)):
        status = "supported" if "supported_in_both" in support_row.external_identity_support else "not supported"
        external_axis.text(1.02, 0.84 - display_index * 0.12, f"{support_row.state}: {status}", transform=external_axis.transAxes, fontsize=5.3, color=COLORS["teal"] if status == "supported" else COLORS["red"])

    compatible = external_candidates.sort_values("minimum_zscore").tail(12).reset_index(drop=True)
    compatible["label"] = compatible["group_members"].str.replace("|", "; ", regex=False).str.slice(0, 38)
    compatible.loc[compatible["group_members"].str.len().gt(38), "label"] += "…"
    positions = np.arange(len(compatible))
    highlight = compatible["group_members"].str.contains("PLAU", regex=False) & compatible["group_members"].str.contains("SERPINE1", regex=False)
    lollipop_axis.hlines(positions, 0, compatible["minimum_zscore"], color="#D5DADF", linewidth=0.8)
    lollipop_axis.scatter(compatible["minimum_zscore"], positions, s=38, c=np.where(highlight, COLORS["vm"], COLORS["teal"]), edgecolor="white", linewidth=0.4, zorder=3)
    lollipop_axis.set_yticks(positions)
    lollipop_axis.set_yticklabels(compatible["label"], fontsize=6.4)
    lollipop_axis.set_xlabel("Minimum cross-cohort receiver-side activity z score")
    lollipop_axis.set_title("External source-RNA-compatible receiver-side associations", loc="left", pad=3)
    lollipop_axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    lollipop_axis.text(0.01, -0.18, "Purple: PLAU/SERPINE1-containing group. Compatibility denotes state-restricted source RNA only; it does not establish secretion, receptor binding, or cell-cell communication.", transform=lollipop_axis.transAxes, fontsize=6.2, va="top")
    panel_label(myeloid_axis, "A")
    panel_label(coverage_axis, "B")
    panel_label(funnel_axis, "C")
    panel_label(external_axis, "D")
    panel_label(lollipop_axis, "E", x=-0.045, y=1.08)
    figure.suptitle("Fig. 3 | Constrained microenvironmental associations linked to the VM-related epithelial programme", x=0.01, y=0.995, ha="left", fontsize=11, fontweight="bold")
    figure.text(0.5, 0.008, "All microenvironmental readouts are association-level and hypothesis-generating. No spatial iCCA colocalisation, molecular secretion or directional communication is claimed.", ha="center", fontsize=6.7)
    myeloid_data.to_csv(SOURCE / "Fig3_corrected_GSE189903_myeloid_embedding.tsv.gz", sep="\t", index=False, compression="gzip")
    eligible_states.to_csv(SOURCE / "Fig3_corrected_myeloid_patient_coverage.tsv", sep="\t", index=False)
    pd.DataFrame({"step": funnel_labels, "n_groups_q80": funnel_values}).to_csv(SOURCE / "Fig3_corrected_secact_funnel_q80.tsv", sep="\t", index=False)
    constrained.to_csv(SECACT / "Table_R3_13_corrected_q80_source_and_external_constraint.tsv", sep="\t", index=False)
    external_data.to_csv(SOURCE / "Fig3_corrected_GSE201425_reference_transfer_embedding.tsv.gz", sep="\t", index=False, compression="gzip")
    compatible.to_csv(SOURCE / "Fig3_corrected_displayed_receiver_associations.tsv", sep="\t", index=False)
    save_figure(figure, "Fig3_corrected_receiver_source_context")


def main() -> None:
    setup_style()
    figure_1()
    figure_2()
    figure_3()


if __name__ == "__main__":
    main()
