"""Rebuild manuscript Figures 2–4 from frozen D-drive analysis results."""
from __future__ import annotations

import gzip
import math
import os
import sys
from pathlib import Path

PROJECT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
PACKAGE_CACHE = PROJECT / "cache" / "python_packages"
sys.path.insert(0, str(PACKAGE_CACHE))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch

RESULTS = PROJECT / "results"
FIGURES = PROJECT / "figures" / "manuscript_v1"
SOURCE_DATA = FIGURES / "source_data"

COLORS = {
    "dark": "#263238",
    "grey": "#9AA5AE",
    "light_grey": "#E9ECEF",
    "vm": "#CC79A7",
    "blue": "#0072B2",
    "teal": "#009E73",
    "orange": "#E69F00",
    "protein": "#D55E00",
    "red": "#A61B1B",
    "pale_blue": "#EAF3F8",
    "pale_pink": "#F9EEF5",
    "pale_green": "#EDF7F2",
    "pale_orange": "#FFF4E5",
}


def setup_style():
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


def save_figure(figure, basename: str):
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURES / f"{basename}.pdf", bbox_inches="tight")
    figure.savefig(FIGURES / f"{basename}.png", dpi=600, bbox_inches="tight")
    plt.close(figure)


def panel_label(axis, label: str, x: float = -0.12, y: float = 1.12):
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


def read_score_file(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def load_embedding(score_data: pd.DataFrame, h5ad_path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(h5ad_path, backed="r")
    coordinates = np.asarray(adata.obsm["X_umap"])
    embedding = pd.DataFrame(
        {"cell_id": adata.obs_names.astype(str), "UMAP_1": coordinates[:, 0], "UMAP_2": coordinates[:, 1]}
    )
    adata.file.close()
    merged = score_data.merge(embedding, on="cell_id", how="left", validate="one_to_one")
    if merged[["UMAP_1", "UMAP_2"]].isna().any().any():
        raise ValueError(f"Embedding lookup failed for {h5ad_path}")
    return merged


def cohort_summary(score_data: pd.DataFrame, cohort: str) -> dict:
    patient_count = score_data["patient_id"].nunique()
    malignant_cells = len(score_data)
    triple_count = int(score_data["VM_triple_consensus_candidate"].astype(bool).sum())
    vm_high_count = int(score_data["VM_high_candidate"].astype(bool).sum())
    return {
        "cohort": cohort,
        "patients": patient_count,
        "operational_malignant_epithelial_cells": malignant_cells,
        "triple_consensus_cells": triple_count,
        "triple_consensus_fraction": triple_count / malignant_cells if malignant_cells else np.nan,
        "VM_high_candidate_cells_after_endothelial_control": vm_high_count,
    }


def plot_embedding(axis, data: pd.DataFrame, title: str):
    background = data.loc[~data["VM_high_candidate"].astype(bool)]
    highlighted = data.loc[data["VM_high_candidate"].astype(bool)]
    axis.scatter(background["UMAP_1"], background["UMAP_2"], s=1.8, c="#C9D0D6", alpha=0.45, linewidths=0)
    axis.scatter(highlighted["UMAP_1"], highlighted["UMAP_2"], s=3.8, c=COLORS["vm"], alpha=0.78, linewidths=0)
    axis.set_title(title, pad=3)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlabel("UMAP 1", labelpad=1)
    axis.set_ylabel("UMAP 2", labelpad=1)


def plot_score_concordance(axis, data: pd.DataFrame, title: str):
    background = data.loc[~data["VM_high_candidate"].astype(bool)]
    highlighted = data.loc[data["VM_high_candidate"].astype(bool)]
    axis.scatter(background["UCell_VM_anchor"], background["AUCell_VM_anchor"], s=2.4, c="#BFC8CF", alpha=0.35, linewidths=0)
    axis.scatter(highlighted["UCell_VM_anchor"], highlighted["AUCell_VM_anchor"], s=4.2, c=COLORS["vm"], alpha=0.75, linewidths=0)
    axis.set_title(title, pad=3)
    axis.set_xlabel("UCell VM-anchor score")
    axis.set_ylabel("AUCell VM-anchor score")
    axis.grid(color=COLORS["light_grey"], linewidth=0.5)


def plot_patient_coverage(axis, data: pd.DataFrame, cohort: str):
    coverage = data.groupby("patient_id", as_index=False).agg(
        operational_malignant_epithelial_cells=("cell_id", "size"),
        triple_consensus_cells=("VM_triple_consensus_candidate", "sum"),
        VM_high_candidate_cells=("VM_high_candidate", "sum"),
    )
    coverage["triple_consensus_fraction"] = coverage["triple_consensus_cells"] / coverage["operational_malignant_epithelial_cells"]
    coverage["cohort"] = cohort
    coverage = coverage.sort_values("triple_consensus_fraction", ascending=True).reset_index(drop=True)
    positions = np.arange(len(coverage))
    scaled_sizes = 14 + 75 * np.sqrt(coverage["operational_malignant_epithelial_cells"] / coverage["operational_malignant_epithelial_cells"].max())
    axis.scatter(coverage["triple_consensus_fraction"] * 100, positions, s=scaled_sizes, c=COLORS["vm"], edgecolor="white", linewidth=0.35, alpha=0.85)
    axis.set_yticks(positions)
    axis.set_yticklabels(coverage["patient_id"].astype(str).str.replace("GSE125449_", "", regex=False), fontsize=5.8)
    axis.set_title(cohort, pad=3)
    axis.set_xlabel("Triple-consensus cells (%)")
    axis.set_xlim(left=0)
    axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    return coverage


def figure_2():
    cohort_inputs = {
        "GSE138709": {
            "scores": RESULTS / "15_gse138709_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
            "h5ad": RESULTS / "11_gse138709_strict_qc_map" / "gse138709_strict_qc_map.h5ad",
        },
        "GSE181878": {
            "scores": RESULTS / "29_gse181878_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
            "h5ad": RESULTS / "26_gse181878_malignant_epithelial_definition" / "gse181878_operational_malignant_epithelium.h5ad",
        },
        "GSE125449": {
            "scores": RESULTS / "46_gse125449_icca_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
            "h5ad": RESULTS / "43_gse125449_icca_malignant_epithelial_definition" / "gse125449_icca_operational_malignant_epithelium.h5ad",
        },
    }
    embedded_data = {}
    summary_rows = []
    for cohort, input_paths in cohort_inputs.items():
        score_data = read_score_file(input_paths["scores"])
        embedded_data[cohort] = load_embedding(score_data, input_paths["h5ad"])
        summary_rows.append(cohort_summary(score_data, cohort))
    cohort_summary_table = pd.DataFrame(summary_rows)
    patient_coverages = []
    figure = plt.figure(figsize=(11.0, 9.4))
    grid = figure.add_gridspec(4, 3, height_ratios=[0.58, 1.08, 1.06, 0.98], hspace=0.70, wspace=0.42)
    flow_axis = figure.add_subplot(grid[0, :])
    flow_axis.axis("off")
    flow_axis.set_xlim(0, 1)
    flow_axis.set_ylim(0, 1)
    flow_axis.text(0.0, 0.98, "A", fontsize=10, fontweight="bold", va="top")
    flow_axis.text(0.07, 0.98, "Cohort-aware operational definition and score rule", fontsize=8.5, fontweight="bold", va="top")
    panel_positions = [0.07, 0.38, 0.69]
    for position, row in zip(panel_positions, cohort_summary_table.itertuples(index=False)):
        flow_axis.add_patch(FancyBboxPatch((position, 0.18), 0.24, 0.56, boxstyle="round,pad=0.018,rounding_size=0.03", facecolor=COLORS["pale_pink"], edgecolor="#D7AEC6", linewidth=0.7))
        flow_axis.text(position + 0.12, 0.63, row.cohort, ha="center", fontsize=8, fontweight="bold")
        patient_label = "patient" if row.patients == 1 else "patients"
        flow_axis.text(position + 0.12, 0.45, f"{row.patients} {patient_label} | {row.operational_malignant_epithelial_cells:,} cells", ha="center", fontsize=6.6)
        flow_axis.text(position + 0.12, 0.31, f"{row.triple_consensus_cells:,} triple-consensus ({row.triple_consensus_fraction * 100:.1f}%)", ha="center", fontsize=6.4, color=COLORS["dark"])
    flow_axis.text(0.50, 0.01, "Operational malignant epithelium → UCell + AUCell + module-score consensus → endothelial-control exclusion", ha="center", fontsize=6.9, color=COLORS["dark"])
    embedding_axes = [figure.add_subplot(grid[1, column]) for column in range(3)]
    concordance_axes = [figure.add_subplot(grid[2, column]) for column in range(3)]
    coverage_axes = [figure.add_subplot(grid[3, column]) for column in range(3)]
    for axis, (cohort, data) in zip(embedding_axes, embedded_data.items()):
        plot_embedding(axis, data, cohort)
    for axis, (cohort, data) in zip(concordance_axes, embedded_data.items()):
        plot_score_concordance(axis, data, cohort)
    for axis, (cohort, data) in zip(coverage_axes, embedded_data.items()):
        patient_coverages.append(plot_patient_coverage(axis, data, cohort))
    panel_label(embedding_axes[0], "B")
    panel_label(concordance_axes[0], "C")
    panel_label(coverage_axes[0], "D")
    figure.suptitle("Fig. 2 | Cross-cohort definition and patient-level coverage of a VM-associated epithelial programme", x=0.01, y=0.995, ha="left", fontsize=11, fontweight="bold")
    figure.text(0.5, 0.012, "Purple cells meet all three prespecified transcriptomic score thresholds and pass the endothelial-control exclusion. They denote a transcriptional state, not morphologic vasculogenic mimicry.", ha="center", fontsize=6.8)
    cohort_summary_table.to_csv(SOURCE_DATA / "Fig2_reconstructed_cohort_summary.tsv", sep="\t", index=False)
    pd.concat(patient_coverages, ignore_index=True).to_csv(SOURCE_DATA / "Fig2_reconstructed_patient_coverage.tsv", sep="\t", index=False)
    combined_embedding = pd.concat([data.assign(cohort=cohort) for cohort, data in embedded_data.items()], ignore_index=True)
    combined_embedding.to_csv(SOURCE_DATA / "Fig2_reconstructed_cell_scores_and_embeddings.tsv.gz", sep="\t", index=False, compression="gzip")
    save_figure(figure, "Fig2_patient_level_VMlike_coverage")


def categorical_embedding(axis, data: pd.DataFrame, state_column: str, title: str, palette: dict):
    categories = [category for category in palette if category in set(data[state_column].astype(str))]
    for category in categories:
        selected = data.loc[data[state_column].astype(str).eq(category)]
        axis.scatter(selected["UMAP_1"], selected["UMAP_2"], s=3.0, color=palette[category], alpha=0.72, linewidths=0, label=category)
    axis.set_title(title, pad=3)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlabel("UMAP 1", labelpad=1)
    axis.set_ylabel("UMAP 2", labelpad=1)


def h5ad_embedding_with_obs(h5ad_path: Path, required_columns: list[str]) -> pd.DataFrame:
    adata = ad.read_h5ad(h5ad_path, backed="r")
    coordinates = np.asarray(adata.obsm["X_umap"])
    selected_obs = adata.obs.loc[:, required_columns].copy()
    selected_obs["cell_id"] = selected_obs.index.astype(str)
    selected_obs["UMAP_1"] = coordinates[:, 0]
    selected_obs["UMAP_2"] = coordinates[:, 1]
    adata.file.close()
    return selected_obs.reset_index(drop=True)


def figure_3():
    myeloid_path = RESULTS / "61_myeloid_stroma_patient_coverage_audit" / "gse189903_myeloid_fine_state_audited.h5ad"
    external_path = RESULTS / "76_gse201425_independent_myeloid_identity_audit" / "gse201425_icca_myeloid_reference_mapped.h5ad"
    myeloid_data = h5ad_embedding_with_obs(myeloid_path, ["fine_state", "case_prefix"])
    external_data = h5ad_embedding_with_obs(external_path, ["reference_state_call", "patient_id"])
    myeloid_coverage = pd.read_csv(RESULTS / "61_myeloid_stroma_patient_coverage_audit" / "myeloid_fine_state_patient_coverage.tsv", sep="\t")
    external_support = pd.read_csv(RESULTS / "76_gse201425_independent_myeloid_identity_audit" / "external_identity_support_summary.tsv", sep="\t")
    status_counts = pd.read_csv(RESULTS / "77_gse201425_secact_patient_state_pseudobulk" / "external_pseudobulk_group_status_counts.tsv", sep="\t")
    candidate_groups = pd.read_csv(RESULTS / "63_secact_cross_cohort_receiver_candidates" / "secact_receiver_activity_nonredundant_groups_cross_cohort.tsv", sep="\t")
    source_restricted = pd.read_csv(RESULTS / "66_secact_state_restricted_myeloid_prioritization" / "secact_receiver_groups_with_state_restricted_myeloid_members.tsv", sep="\t")
    multilayer = pd.read_csv(RESULTS / "78_secact_multilayer_evidence_matrix" / "secact_30_multilayer_evidence_matrix.tsv", sep="\t")
    candidate_count = len(candidate_groups)
    positive_count = int(candidate_groups["cross_cohort_signature_candidate"].astype(bool).sum())
    source_count = int(source_restricted["has_state_restricted_myeloid_source_member"].astype(bool).sum())
    external_count = int(status_counts.loc[status_counts["status"].eq("has_externally_consistent_source_RNA_member"), "n_groups"].iloc[0])
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
    legend_handles, legend_labels = myeloid_axis.get_legend_handles_labels()
    eligible_state_labels = set(eligible_states["fine_state"].astype(str))
    eligible_handles = [handle for handle, label in zip(legend_handles, legend_labels) if label in eligible_state_labels]
    eligible_labels = [label.replace("_", " ") for label in legend_labels if label in eligible_state_labels]
    myeloid_axis.legend(eligible_handles, eligible_labels, frameon=False, markerscale=2.0, fontsize=4.7, loc="upper center", bbox_to_anchor=(0.50, -0.19), ncol=2)
    states_for_heatmap = eligible_states["fine_state"].astype(str).tolist()
    heatmap_values = eligible_states.loc[:, ["cells_1C", "cells_2C", "cells_3C"]].to_numpy(dtype=float)
    image = coverage_axis.imshow(heatmap_values, cmap="Purples", aspect="auto")
    coverage_axis.set_xticks(range(3))
    coverage_axis.set_xticklabels(["Patient 1C", "Patient 2C", "Patient 3C"])
    coverage_axis.set_yticks(range(len(states_for_heatmap)))
    coverage_axis.set_yticklabels([state.replace("_", " ") for state in states_for_heatmap], fontsize=5.8)
    coverage_axis.set_title("Patient coverage of eligible reference states", pad=3)
    for row_index in range(heatmap_values.shape[0]):
        for column_index in range(heatmap_values.shape[1]):
            coverage_axis.text(column_index, row_index, f"{int(heatmap_values[row_index, column_index])}", ha="center", va="center", fontsize=6.2, color="white" if heatmap_values[row_index, column_index] > np.nanmedian(heatmap_values) else COLORS["dark"])
    coverage_axis.figure.colorbar(image, ax=coverage_axis, fraction=0.045, pad=0.02, label="Cells")
    funnel_labels = ["Nonredundant SecAct groups", "Positive in both scRNA cohorts", "State-restricted myeloid-source member", "External patient×state RNA-compatible"]
    funnel_values = [candidate_count, positive_count, source_count, external_count]
    funnel_colors = ["#B7C3CE", COLORS["blue"], COLORS["teal"], COLORS["vm"]]
    funnel_axis.barh(np.arange(len(funnel_values))[::-1], funnel_values, color=funnel_colors, height=0.62)
    funnel_axis.set_yticks(np.arange(len(funnel_values))[::-1])
    funnel_axis.set_yticklabels(funnel_labels, fontsize=6.1)
    funnel_axis.set_xlabel("Number of groups")
    funnel_axis.set_title("Predefined evidence-constraining workflow", pad=3)
    for display_index, value in enumerate(funnel_values[::-1]):
        funnel_axis.text(value + max(funnel_values) * 0.012, display_index, str(value), va="center", fontsize=6.8, fontweight="bold")
    funnel_axis.set_xlim(0, max(funnel_values) * 1.14)
    funnel_axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    categorical_embedding(external_axis, external_data, "reference_state_call", "GSE201425 independent reference transfer", external_palette)
    support_text = external_support.assign(state=external_support["state"].str.replace("_", " ", regex=False)).loc[:, ["state", "external_identity_support"]]
    external_axis.text(1.02, 0.98, "External identity support", transform=external_axis.transAxes, va="top", fontsize=6.2, fontweight="bold")
    for display_index, support_row in enumerate(support_text.itertuples(index=False)):
        status = "supported" if "supported_in_both" in support_row.external_identity_support else "not supported"
        external_axis.text(1.02, 0.84 - display_index * 0.12, f"{support_row.state}: {status}", transform=external_axis.transAxes, fontsize=5.3, color=COLORS["teal"] if status == "supported" else COLORS["red"])
    compatible = multilayer.loc[multilayer["external_pseudobulk_group_status"].eq("has_externally_consistent_source_RNA_member")].copy()
    compatible = compatible.sort_values("minimum_zscore", ascending=True).tail(12).reset_index(drop=True)
    compatible["label"] = compatible["group_members"].str.replace("|", "; ", regex=False).str.slice(0, 38)
    compatible.loc[compatible["group_members"].str.len().gt(38), "label"] = compatible.loc[compatible["group_members"].str.len().gt(38), "label"] + "…"
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
    figure.suptitle("Fig. 3 | Constrained microenvironmental associations linked to the VM-like epithelial state", x=0.01, y=0.995, ha="left", fontsize=11, fontweight="bold")
    figure.text(0.5, 0.008, "All microenvironmental readouts are association-level and hypothesis-generating. No spatial iCCA colocalization, molecular secretion, or directional communication is claimed.", ha="center", fontsize=6.7)
    myeloid_data.to_csv(SOURCE_DATA / "Fig3_reconstructed_GSE189903_myeloid_embedding.tsv.gz", sep="\t", index=False, compression="gzip")
    eligible_states.to_csv(SOURCE_DATA / "Fig3_reconstructed_myeloid_patient_coverage.tsv", sep="\t", index=False)
    pd.DataFrame({"step": funnel_labels, "n_groups": funnel_values}).to_csv(SOURCE_DATA / "Fig3_reconstructed_secact_funnel.tsv", sep="\t", index=False)
    external_data.to_csv(SOURCE_DATA / "Fig3_reconstructed_GSE201425_reference_transfer_embedding.tsv.gz", sep="\t", index=False, compression="gzip")
    compatible.to_csv(SOURCE_DATA / "Fig3_reconstructed_displayed_receiver_associations.tsv", sep="\t", index=False)
    save_figure(figure, "Fig3_constrained_receiver_associations")


def forest_axis(axis, selected: pd.DataFrame):
    displayed = selected.sort_values("fixed_beta", ascending=True).reset_index(drop=True)
    positions = np.arange(len(displayed))
    cohort_columns = {"GSE32225_beta": (-0.20, COLORS["blue"], "GSE32225"), "GSE107943_beta": (0.0, COLORS["orange"], "GSE107943"), "GSE26566_beta": (0.20, COLORS["teal"], "GSE26566")}
    for column, (offset, color, label) in cohort_columns.items():
        axis.scatter(displayed[column], positions + offset, s=22, color=color, label=label, zorder=3)
    x_errors = [displayed["fixed_beta"] - displayed["fixed_ci_low"], displayed["fixed_ci_high"] - displayed["fixed_beta"]]
    axis.errorbar(displayed["fixed_beta"], positions, xerr=x_errors, fmt="s", color=COLORS["dark"], markersize=3.8, capsize=2, label="Three-cohort fixed effect", zorder=4)
    axis.axvline(0, color="#808A93", linewidth=0.7)
    axis.set_yticks(positions)
    axis.set_yticklabels(displayed["external_source_RNA_members"].str.replace(";", "; ", regex=False), fontsize=6.4)
    axis.set_xlabel("Association with bulk VM-like score (β)")
    axis.set_title("Three-cohort directional consistency", loc="left", pad=3)
    axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    axis.legend(frameon=False, ncol=2, fontsize=5.9, loc="lower right")
    return displayed


def figure_4():
    meta = pd.read_csv(RESULTS / "87_three_bulk_cohort_meta_and_triage" / "three_bulk_cohort_effect_summary_and_meta.tsv", sep="\t")
    selected = meta.loc[meta["all_three_direction_consistent"].astype(bool)].copy()
    discovery = pd.read_csv(RESULTS / "83_bulk_stromal_endothelial_adjustment" / "vm_associations_with_myeloid_fibroblast_endothelial_controls.tsv", sep="\t")
    validation = pd.read_csv(RESULTS / "84_gse107943_direction_and_clinical_sensitivity" / "seven_group_directional_replication.tsv", sep="\t")
    third = pd.read_csv(RESULTS / "86_gse26566_third_bulk_directional_replication" / "seven_group_third_bulk_directional_replication.tsv", sep="\t")
    selected_groups = selected["signature_group"].tolist()
    discovery = discovery.loc[discovery["signature_group"].isin(selected_groups)].copy()
    validation = validation.loc[validation["signature_group"].isin(selected_groups)].copy()
    third = third.loc[third["signature_group"].isin(selected_groups)].copy()
    sensitivity = selected.loc[:, ["signature_group", "external_source_RNA_members"]].copy()
    for cohort_name, source, adjusted_column, mmp2_column in [
        ("GSE32225", discovery, "four_composition_controls_beta", "without_MMP2_four_controls_beta"),
        ("GSE107943", validation, "gse107943_four_controls_beta", "gse107943_without_MMP2_four_controls_beta"),
        ("GSE26566", third, "GSE26566_four_controls_beta", "GSE26566_without_MMP2_four_controls_beta"),
    ]:
        subset = source.loc[:, ["signature_group", adjusted_column, mmp2_column]].rename(columns={adjusted_column: f"{cohort_name}_four_proxy_beta", mmp2_column: f"{cohort_name}_without_MMP2_beta"})
        sensitivity = sensitivity.merge(subset, on="signature_group", how="left", validate="one_to_one")
    sensitivity = sensitivity.set_index("signature_group").loc[selected_groups].reset_index()
    heatmap_columns = ["GSE32225_four_proxy_beta", "GSE32225_without_MMP2_beta", "GSE107943_four_proxy_beta", "GSE107943_without_MMP2_beta", "GSE26566_four_proxy_beta", "GSE26566_without_MMP2_beta"]
    heatmap_values = sensitivity.loc[:, heatmap_columns].to_numpy(dtype=float)
    figure = plt.figure(figsize=(11.0, 8.5))
    grid = figure.add_gridspec(2, 2, height_ratios=[1.08, 0.92], hspace=0.58, wspace=0.50)
    forest = figure.add_subplot(grid[0, :])
    heatmap = figure.add_subplot(grid[1, 0])
    meta_axis = figure.add_subplot(grid[1, 1])
    displayed = forest_axis(forest, selected)
    image = heatmap.imshow(heatmap_values, cmap="RdBu_r", vmin=-0.75, vmax=0.75, aspect="auto")
    heatmap.set_xticks(range(len(heatmap_columns)))
    heatmap.set_xticklabels(["32225\n4 proxies", "32225\n−MMP2", "107943\n4 proxies", "107943\n−MMP2", "26566\n4 proxies", "26566\n−MMP2"], fontsize=5.7)
    heatmap.set_yticks(range(len(sensitivity)))
    heatmap.set_yticklabels(sensitivity["external_source_RNA_members"].str.replace(";", "; ", regex=False), fontsize=5.8)
    heatmap.set_title("Composition and anchor sensitivity", loc="left", pad=3)
    for row_index in range(heatmap_values.shape[0]):
        for column_index in range(heatmap_values.shape[1]):
            value = heatmap_values[row_index, column_index]
            heatmap.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=5.2, color="white" if abs(value) > 0.40 else COLORS["dark"])
    heatmap.figure.colorbar(image, ax=heatmap, fraction=0.046, pad=0.03, label="β")
    meta_display = selected.sort_values("external_only_fixed_beta", ascending=True).reset_index(drop=True)
    positions = np.arange(len(meta_display))
    all_errors = [meta_display["fixed_beta"] - meta_display["fixed_ci_low"], meta_display["fixed_ci_high"] - meta_display["fixed_beta"]]
    external_errors = [meta_display["external_only_fixed_beta"] - meta_display["external_only_fixed_ci_low"], meta_display["external_only_fixed_ci_high"] - meta_display["external_only_fixed_beta"]]
    meta_axis.errorbar(meta_display["fixed_beta"], positions + 0.14, xerr=all_errors, fmt="s", color=COLORS["dark"], markersize=3.7, capsize=2, label="Discovery-inclusive fixed effect")
    meta_axis.errorbar(meta_display["external_only_fixed_beta"], positions - 0.14, xerr=external_errors, fmt="o", color=COLORS["vm"], markersize=3.7, capsize=2, label="External-only fixed effect")
    meta_axis.axvline(0, color="#808A93", linewidth=0.7)
    meta_axis.set_yticks(positions)
    meta_axis.set_yticklabels(meta_display["external_source_RNA_members"].str.replace(";", "; ", regex=False), fontsize=5.8)
    meta_axis.set_xlabel("Meta-analytic β")
    meta_axis.set_title("Discovery-inclusive versus external-only synthesis", loc="left", pad=3)
    meta_axis.legend(frameon=False, fontsize=5.6, loc="lower right")
    meta_axis.grid(axis="x", color=COLORS["light_grey"], linewidth=0.5)
    panel_label(forest, "A", x=-0.055)
    panel_label(heatmap, "B", x=-0.14)
    panel_label(meta_axis, "C", x=-0.14)
    figure.suptitle("Fig. 4 | Directional bulk consistency and composition-aware sensitivity across iCCA cohorts", x=0.01, y=0.995, ha="left", fontsize=11, fontweight="bold")
    figure.text(0.5, 0.012, "All models are patient-sample-level bulk associations adjusted with expression proxies for myeloid, fibroblast and endothelial composition; they do not identify a cellular source or establish causality.", ha="center", fontsize=6.7)
    displayed.to_csv(SOURCE_DATA / "Fig4_reconstructed_three_cohort_directional_replication.tsv", sep="\t", index=False)
    sensitivity.to_csv(SOURCE_DATA / "Fig4_reconstructed_composition_and_MMP2_sensitivity.tsv", sep="\t", index=False)
    meta_display.to_csv(SOURCE_DATA / "Fig4_reconstructed_discovery_and_external_meta.tsv", sep="\t", index=False)
    save_figure(figure, "Fig4_three_bulk_directional_replication")


def main():
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    setup_style()
    figure_2()
    figure_3()
    figure_4()


if __name__ == "__main__":
    main()
