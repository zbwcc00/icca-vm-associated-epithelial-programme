from pathlib import Path
import gzip
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
RESULTS = PROJECT / "results"
OUTPUT = PROJECT / "figures" / "manuscript_v1"
SOURCE_DATA = OUTPUT / "source_data"
ANCHORS = ["CDH5", "EPHA2", "MCAM", "MMP2", "MMP9", "LAMC2"]
COLORS = {
    "vm": "#CC79A7",
    "rna": "#0072B2",
    "protein": "#D55E00",
    "teal": "#009E73",
    "orange": "#E69F00",
    "grey": "#6C757D",
    "light": "#E9ECEF",
    "dark": "#343A40",
}


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save_figure(figure, filename):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / f"{filename}.pdf", bbox_inches="tight")
    figure.savefig(OUTPUT / f"{filename}.png", dpi=600, bbox_inches="tight")
    plt.close(figure)


def panel_label(axis, label, x=-0.08, y=1.11):
    axis.text(x, y, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top", bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.15})


def draw_box(axis, xy, width, height, text, facecolor, edgecolor=COLORS["dark"]):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.015,rounding_size=0.02", linewidth=0.8, facecolor=facecolor, edgecolor=edgecolor)
    axis.add_patch(box)
    axis.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=7, wrap=True)


def figure_1():
    figure, axis = plt.subplots(figsize=(7.7, 5.1))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.add_patch(FancyBboxPatch((0.015, 0.015), 0.97, 0.97, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=0.6, facecolor="#FCFCFD", edgecolor="#D9DEE5"))
    axis.text(0.04, 0.945, "Fig. 1 | Evidence architecture of a VM-associated epithelial programme", fontsize=10.5, fontweight="bold", va="center")
    axis.text(0.04, 0.912, "Public iCCA cohorts are integrated in a patient-aware framework; arrows show evidence flow, not causality.", fontsize=7, color=COLORS["grey"], va="center")
    axis.plot([0.04, 0.96], [0.885, 0.885], color="#D9DEE5", linewidth=0.8)

    def card(x, y, width, height, title, body, accent, dashed=False):
        line_style = "--" if dashed else "-"
        axis.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.018", linewidth=0.9, linestyle=line_style, facecolor="#FFFFFF", edgecolor=accent, zorder=2))
        axis.add_patch(FancyBboxPatch((x, y + height - 0.043), width, 0.043, boxstyle="round,pad=0.012,rounding_size=0.018", linewidth=0, facecolor=accent, edgecolor=accent, zorder=3))
        axis.add_patch(Circle((x + 0.022, y + height - 0.0215), 0.008, facecolor="#FFFFFF", edgecolor="none", zorder=4))
        axis.text(x + 0.04, y + height - 0.0215, title, fontsize=6.2, fontweight="bold", color="#FFFFFF", va="center", zorder=4)
        axis.text(x + width / 2, y + (height - 0.043) / 2, body, fontsize=6.1, ha="center", va="center", linespacing=1.35, zorder=4)

    card(0.045, 0.62, 0.245, 0.19, "01  SINGLE-CELL MAP", "3 iCCA scRNA-seq cohorts\n15 patients; malignant epithelium\ntriple-score coverage", COLORS["rna"])
    card(0.71, 0.62, 0.245, 0.19, "02  BULK CONSISTENCY", "3 independent bulk cohorts\n6 directionally concordant groups\ncomposition-aware sensitivity", COLORS["teal"])
    card(0.045, 0.305, 0.245, 0.17, "03  HYPOTHESIS CONTEXT", "Receiver activity plus\nexternal source-RNA compatibility\nC1QC-TAM candidate context", COLORS["orange"], dashed=True)
    card(0.71, 0.305, 0.245, 0.17, "04  MATCHED OMICS", "OEP001105 iCCA tumours\nRNA–protein: n=208; ρ=0.709\nRNA–WES/CNV: n=249", COLORS["protein"])

    ellipse = Ellipse((0.5, 0.57), 0.35, 0.25, facecolor="#F9D8E8", edgecolor=COLORS["vm"], linewidth=1.4, zorder=3)
    axis.add_patch(ellipse)
    axis.text(0.5, 0.627, "VM-associated malignant\nepithelial programme", ha="center", va="center", fontsize=8.5, fontweight="bold", color="#5D254C", zorder=4)
    axis.text(0.5, 0.552, "six prespecified anchors", ha="center", va="center", fontsize=6.5, color="#5D254C", zorder=4)
    for x, label in zip([0.425, 0.5, 0.575], ["CDH5", "EPHA2", "MCAM"]):
        axis.add_patch(FancyBboxPatch((x - 0.031, 0.50), 0.062, 0.027, boxstyle="round,pad=0.004,rounding_size=0.012", linewidth=0.3, facecolor="#FFFFFF", edgecolor="#D89ABB", zorder=4))
        axis.text(x, 0.5135, label, ha="center", va="center", fontsize=5.6, color="#5D254C", zorder=5)
    for x, label in zip([0.425, 0.5, 0.575], ["MMP2", "MMP9", "LAMC2"]):
        axis.add_patch(FancyBboxPatch((x - 0.031, 0.465), 0.062, 0.027, boxstyle="round,pad=0.004,rounding_size=0.012", linewidth=0.3, facecolor="#FFFFFF", edgecolor="#D89ABB", zorder=4))
        axis.text(x, 0.4785, label, ha="center", va="center", fontsize=5.6, color="#5D254C", zorder=5)

    def flow(start, end, color, linestyle="-"):
        axis.add_patch(FancyArrowPatch(start, end, connectionstyle="arc3,rad=0.12", arrowstyle="Simple,tail_width=0.35,head_width=4,head_length=5", linewidth=0.5, linestyle=linestyle, color=color, alpha=0.85, zorder=1))

    flow((0.29, 0.69), (0.355, 0.64), COLORS["rna"])
    flow((0.71, 0.69), (0.645, 0.64), COLORS["teal"])
    flow((0.29, 0.385), (0.36, 0.48), COLORS["orange"], "--")
    flow((0.71, 0.385), (0.64, 0.48), COLORS["protein"])
    axis.text(0.5, 0.395, "Cross-scale integration", ha="center", va="center", fontsize=7, color=COLORS["grey"])

    axis.add_patch(FancyBboxPatch((0.045, 0.095), 0.91, 0.12, boxstyle="round,pad=0.014,rounding_size=0.018", linewidth=0, facecolor="#EEF1F4", zorder=2))
    axis.text(0.065, 0.172, "FRAMEWORK-LEVEL CONCLUSION", fontsize=7, fontweight="bold", color=COLORS["dark"], va="center")
    axis.text(0.065, 0.145, "VM-associated epithelial programme observed across heterogeneous iCCA cohorts.", fontsize=6.55, color=COLORS["dark"], va="center")
    axis.text(0.065, 0.118, "Within-project RNA–protein concordance; myeloid and genomic association-level context.", fontsize=6.55, color=COLORS["dark"], va="center")
    axis.text(0.5, 0.045, "Boundary: no vasculogenic-mimicry channels/perfusion, secretion/communication, or phosphorylation activation;\nno treatment resistance, independent prognosis, or drug vulnerability.", ha="center", va="center", fontsize=5.65, color="#9B1C1C")
    save_figure(figure, "Fig1_study_design_and_claim_boundaries")


def read_score_file(path):
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def figure_2():
    input_paths = {
        "GSE138709": RESULTS / "15_gse138709_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
        "GSE181878": RESULTS / "29_gse181878_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
        "GSE125449": RESULTS / "46_gse125449_icca_vm_triple_scoring" / "malignant_epithelial_vm_triple_scores.tsv.gz",
    }
    summaries = []
    for cohort, path in input_paths.items():
        cells = read_score_file(path)
        patient = cells.groupby("patient_id", as_index=False).agg(
            malignant_epithelial_cells=("cell_id", "size"),
            vm_triple_consensus_cells=("VM_triple_consensus_candidate", "sum"),
        )
        patient["vm_triple_consensus_fraction"] = patient["vm_triple_consensus_cells"] / patient["malignant_epithelial_cells"]
        patient["cohort"] = cohort
        summaries.append(patient)
    coverage = pd.concat(summaries, ignore_index=True)
    coverage.to_csv(SOURCE_DATA / "Fig2_patient_level_VMlike_coverage.tsv", sep="\t", index=False)
    figure, axes = plt.subplots(1, 3, figsize=(7.1, 4.2), sharex=True)
    for axis, cohort in zip(axes, input_paths):
        data = coverage.loc[coverage["cohort"].eq(cohort)].sort_values("vm_triple_consensus_fraction", ascending=True).reset_index(drop=True)
        y = np.arange(len(data))
        sizes = 16 + 90 * np.sqrt(data["malignant_epithelial_cells"] / data["malignant_epithelial_cells"].max())
        axis.scatter(data["vm_triple_consensus_fraction"] * 100, y, s=sizes, color=COLORS["vm"], alpha=0.8, edgecolor="white", linewidth=0.4)
        axis.set_yticks(y)
        axis.set_yticklabels(data["patient_id"].str.replace("GSE125449_", "", regex=False), fontsize=6)
        axis.set_title(f"{cohort} (n={len(data)} patients)")
        axis.grid(axis="x", color="#E9ECEF", linewidth=0.7)
        axis.set_xlabel("Triple-consensus VM-like cells (%)")
        axis.set_xlim(left=0)
    axes[0].set_ylabel("Patient")
    panel_label(axes[0], "A")
    figure.suptitle("Fig. 2 | Patient-level coverage of VM-like malignant epithelial cells", x=0.01, y=1.02, ha="left", fontsize=10, fontweight="bold")
    figure.text(0.5, -0.04, "Dot size denotes retained malignant epithelial cells; fractions are descriptive and do not establish morphologic VM.", ha="center", fontsize=7)
    figure.tight_layout()
    save_figure(figure, "Fig2_patient_level_VMlike_coverage")


def figure_3():
    source = RESULTS / "78_secact_multilayer_evidence_matrix" / "secact_30_multilayer_evidence_matrix.tsv"
    data = pd.read_csv(source, sep="\t")
    compatible = data.loc[data["external_pseudobulk_group_status"].eq("has_externally_consistent_source_RNA_member")].copy()
    if compatible.empty:
        compatible = data.loc[data["n_members_with_external_state_restricted_source_RNA"].gt(0)].copy()
    selected = compatible.sort_values("minimum_zscore", ascending=False).head(12).copy()
    selected["display"] = selected["group_members"].str.replace("|", "; ", regex=False).str.slice(0, 34)
    selected.loc[selected["group_members"].str.len().gt(34), "display"] = selected.loc[selected["group_members"].str.len().gt(34), "display"] + "…"
    selected.to_csv(SOURCE_DATA / "Fig3_receiver_activity_source_RNA_compatibility.tsv", sep="\t", index=False)
    figure, axis = plt.subplots(figsize=(7.1, 4.5))
    selected = selected.sort_values("minimum_zscore", ascending=True).reset_index(drop=True)
    y = np.arange(len(selected))
    is_plau = selected["group_members"].str.contains("PLAU", regex=False) & selected["group_members"].str.contains("SERPINE1", regex=False)
    colors = np.where(is_plau, COLORS["vm"], COLORS["teal"])
    sizes = 45 + 45 * selected["n_members_with_external_state_restricted_source_RNA"].fillna(0)
    axis.scatter(selected["minimum_zscore"], y, s=sizes, c=colors, edgecolor="white", linewidth=0.5, zorder=3)
    axis.hlines(y, 0, selected["minimum_zscore"], color="#CED4DA", linewidth=0.9, zorder=1)
    axis.set_yticks(y)
    axis.set_yticklabels(selected["display"])
    axis.set_xlabel("Minimum cross-cohort receiver-activity z score")
    axis.set_ylabel("")
    axis.grid(axis="x", color="#E9ECEF", linewidth=0.7)
    axis.set_title("Fig. 3 | Constrained receiver-side associations with external source-RNA compatibility", loc="left", fontweight="bold")
    axis.text(0.01, -0.18, "Purple: PLAU/SERPINE1-containing group. Dot size: number of externally consistent source-RNA members. This is not evidence of secretion, ligand-receptor binding, or cell-cell communication.", transform=axis.transAxes, fontsize=7, va="top", wrap=True)
    save_figure(figure, "Fig3_constrained_receiver_associations")


def figure_4():
    source = RESULTS / "87_three_bulk_cohort_meta_and_triage" / "three_bulk_cohort_effect_summary_and_meta.tsv"
    data = pd.read_csv(source, sep="\t")
    consistent = data["all_three_direction_consistent"].astype(str).str.lower().eq("true")
    selected = data.loc[consistent].copy().sort_values("fixed_beta", ascending=True)
    if len(selected) > 7:
        selected = selected.tail(7)
    selected["display"] = selected["external_source_RNA_members"].fillna(selected["signature_group"])
    selected.to_csv(SOURCE_DATA / "Fig4_three_bulk_directional_replication.tsv", sep="\t", index=False)
    figure, axis = plt.subplots(figsize=(7.1, 4.5))
    y = np.arange(len(selected))
    offsets = {"GSE32225_beta": -0.18, "GSE107943_beta": 0, "GSE26566_beta": 0.18}
    colors = {"GSE32225_beta": COLORS["rna"], "GSE107943_beta": COLORS["orange"], "GSE26566_beta": COLORS["teal"]}
    labels = {"GSE32225_beta": "GSE32225", "GSE107943_beta": "GSE107943", "GSE26566_beta": "GSE26566"}
    for column in offsets:
        axis.scatter(selected[column], y + offsets[column], s=26, color=colors[column], label=labels[column], zorder=3)
    axis.errorbar(selected["fixed_beta"], y, xerr=[selected["fixed_beta"] - selected["fixed_ci_low"], selected["fixed_ci_high"] - selected["fixed_beta"]], fmt="s", color=COLORS["dark"], markersize=4, capsize=2, label="Fixed-effect pooled", zorder=4)
    axis.axvline(0, color=COLORS["grey"], linewidth=0.8)
    axis.set_yticks(y)
    axis.set_yticklabels(selected["display"])
    axis.set_xlabel("Association with bulk VM-like score (regression coefficient)")
    axis.set_ylabel("SecAct group / external source-RNA member")
    axis.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), borderaxespad=0)
    axis.grid(axis="x", color="#E9ECEF", linewidth=0.7)
    axis.set_title("Fig. 4 | Directional replication across three iCCA bulk cohorts", loc="left", fontweight="bold")
    axis.text(0.01, -0.18, "All displayed groups had directionally consistent estimates across the three cohorts. Bulk associations cannot establish cellular source or causality.", transform=axis.transAxes, fontsize=7, va="top", wrap=True)
    save_figure(figure, "Fig4_three_bulk_directional_replication")


def figure_5():
    shared_source = RESULTS / "104_oep001105_multioics_state_validation" / "shared_sample_vm_and_endothelial_scores.tsv"
    anchor_source = RESULTS / "104_oep001105_multioics_state_validation" / "vm_anchor_gene_cross_omic_correlations.tsv"
    clinical_source = RESULTS / "104_oep001105_multioics_state_validation" / "oep001105_vm_clinical_association_tests.tsv"
    shared = pd.read_csv(shared_source, sep="\t")
    shared = shared.dropna(subset=["VM_RNA_log2TPM", "VM_protein_NA"]).copy()
    anchors = pd.read_csv(anchor_source, sep="\t")
    anchors = anchors.loc[(anchors["left_layer"].eq("RNA_log2TPM")) & (anchors["right_layer"].eq("protein_NA")) & anchors["rho"].notna()].copy()
    clinical = pd.read_csv(clinical_source, sep="\t")
    clinical = clinical.loc[clinical["outcome"].eq("VM_RNA_log2TPM_255")].copy()
    shared.to_csv(SOURCE_DATA / "Fig5_matched_RNA_protein_program_scores.tsv", sep="\t", index=False)
    anchors.to_csv(SOURCE_DATA / "Fig5_anchor_RNA_protein_correlations.tsv", sep="\t", index=False)
    clinical.to_csv(SOURCE_DATA / "Fig5_available_clinical_covariate_tests.tsv", sep="\t", index=False)
    rho, p_value = spearmanr(shared["VM_RNA_log2TPM"], shared["VM_protein_NA"])
    figure, axes = plt.subplots(1, 3, figsize=(7.5, 3.2), gridspec_kw={"width_ratios": [1.2, 0.9, 0.9]})
    axes[0].scatter(shared["VM_RNA_log2TPM"], shared["VM_protein_NA"], s=12, alpha=0.65, color=COLORS["vm"], edgecolor="none")
    coefficients = np.polyfit(shared["VM_RNA_log2TPM"], shared["VM_protein_NA"], 1)
    x_values = np.linspace(shared["VM_RNA_log2TPM"].min(), shared["VM_RNA_log2TPM"].max(), 100)
    axes[0].plot(x_values, np.polyval(coefficients, x_values), color=COLORS["dark"], linewidth=1)
    axes[0].set_xlabel("RNA VM-like score")
    axes[0].set_ylabel("Matched total-protein VM-like score")
    axes[0].text(0.02, 1.015, f"n={len(shared)}\nSpearman ρ={rho:.3f}\nP={p_value:.2e}", transform=axes[0].transAxes, va="bottom", fontsize=6.5,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 1.0, "pad": 0.25}, clip_on=False)
    axes[0].grid(color="#E9ECEF", linewidth=0.7)
    anchors = anchors.sort_values("rho")
    axes[1].barh(anchors["gene"], anchors["rho"], color=COLORS["protein"])
    axes[1].set_xlabel("RNA–protein Spearman ρ")
    axes[1].set_xlim(0, 1)
    axes[1].grid(axis="x", color="#E9ECEF", linewidth=0.7)
    clinical["label"] = clinical["test"].str.replace("Spearman: ", "", regex=False).str.replace("Mann-Whitney: ", "", regex=False)
    clinical = clinical.sort_values("fdr_bh_all_prespecified_tests")
    values = -np.log10(clinical["fdr_bh_all_prespecified_tests"].clip(lower=1e-300))
    axes[2].scatter(values, np.arange(len(clinical)), s=28, color=COLORS["grey"])
    axes[2].axvline(-math.log10(0.05), color="#8B0000", linestyle="--", linewidth=0.8)
    axes[2].set_yticks(np.arange(len(clinical)))
    axes[2].set_yticklabels(clinical["label"])
    axes[2].set_xlabel("−log10 BH-FDR")
    axes[2].grid(axis="x", color="#E9ECEF", linewidth=0.7)
    panel_label(axes[0], "A", x=-0.20, y=1.24)
    panel_label(axes[1], "B", x=-0.20, y=1.24)
    panel_label(axes[2], "C", x=-0.20, y=1.24)
    figure.suptitle("Fig. 5 | Orthogonal proteomic support for the VM-like epithelial state", x=0.01, y=1.12, ha="left", fontsize=10, fontweight="bold")
    figure.text(0.5, -0.05, "Phosphopeptide data are limited to EPHA2, MCAM and LAMC2 and are deliberately not combined into a phospho-VM score.", ha="center", fontsize=7)
    figure.tight_layout(w_pad=2.5)
    save_figure(figure, "Fig5_matched_proteomic_concordance")


def functional_mutations(mutations):
    nonfunctional = {
        "synonymous_variant",
        "intron_variant",
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intergenic_variant",
    }
    return mutations.loc[~mutations["Mutation_Type"].isin(nonfunctional)].copy()


def normalize_id(value):
    value = str(value).strip()
    return value[:-2] if value.endswith(".0") else value


def RNA_vm_score(rna_path):
    rna = pd.read_csv(rna_path, sep="\t")
    rna.columns = [normalize_id(column) for column in rna.columns]
    sample_columns = [column for column in rna.columns if column != "GeneSymbol"]
    matrix = rna.set_index("GeneSymbol").apply(pd.to_numeric, errors="coerce")
    matrix.index = matrix.index.astype(str).str.upper()
    matrix = matrix.groupby(level=0).median()
    anchor_matrix = np.log2(matrix.loc[ANCHORS, sample_columns] + 1)
    standardized = anchor_matrix.sub(anchor_matrix.mean(axis=1), axis=0).div(anchor_matrix.std(axis=1, ddof=0), axis=0)
    return standardized.mean(axis=0, skipna=True).rename("VM_RNA_log2TPM")


def figure_6():
    rna_source = PROJECT / "data" / "00_inbox" / "08_OEP001105_proteogenomics" / "01_RNA_processed" / "RNAseq_TPM_255_samples_CPTAC_iCCA.txt"
    snv_source = PROJECT / "data" / "00_inbox" / "08_OEP001105_proteogenomics" / "04_WES_processed" / "WES_SNV_Somatic_mutation_253_samples_CPTAC_iCCA.txt"
    cnv_source = PROJECT / "data" / "00_inbox" / "08_OEP001105_proteogenomics" / "04_WES_processed" / "WES_CNV_GISTIC2_cutoff.txt"
    rna_score = RNA_vm_score(rna_source)
    shared = rna_score.rename_axis("sample_id").reset_index()
    mutations = functional_mutations(pd.read_csv(snv_source, sep="\t", low_memory=False))
    mutations["Sample_ID"] = mutations["Sample_ID"].map(normalize_id)
    mutations["Gene"] = mutations["Gene"].astype(str).str.upper()
    burden = mutations.groupby("Sample_ID").size().rename("nonsynonymous_mutation_burden").reset_index().rename(columns={"Sample_ID": "sample_id"})
    cnv = pd.read_csv(cnv_source, sep="\t", dtype=str)
    cnv = cnv.rename(columns={cnv.columns[0]: "gene_locus"})
    cnv["gene"] = cnv["gene_locus"].str.split("|").str[0].str.upper()
    cnv = cnv.drop_duplicates("gene", keep="first").set_index("gene")
    sample_columns = [column for column in cnv.columns if normalize_id(column) in set(shared["sample_id"])]
    cnv_matrix = cnv.loc[:, sample_columns].apply(pd.to_numeric, errors="coerce")
    cnv_matrix.columns = [normalize_id(column) for column in cnv_matrix.columns]
    cnv_burden = pd.DataFrame({"sample_id": cnv_matrix.columns, "non_diploid_CNV_fraction": cnv_matrix.ne(0).mean(axis=0).values})
    genomic = shared.merge(burden, on="sample_id", how="inner").merge(cnv_burden, on="sample_id", how="inner")
    for gene in ["KRAS", "BAP1"]:
        altered = set(mutations.loc[mutations["Gene"].eq(gene), "Sample_ID"])
        genomic[f"{gene}_mutant"] = genomic["sample_id"].isin(altered)
    for gene in ["MYC", "MCL1"]:
        gain_samples = set(cnv_matrix.columns[cnv_matrix.loc[gene].gt(0)])
        genomic[f"{gene}_gain"] = genomic["sample_id"].isin(gain_samples)
    genomic.to_csv(SOURCE_DATA / "Fig6_matched_genomic_context.tsv", sep="\t", index=False)
    figure = plt.figure(figsize=(7.5, 5.0))
    grid = figure.add_gridspec(2, 4, height_ratios=[1, 1], wspace=0.85, hspace=0.58)
    top_left = figure.add_subplot(grid[0, 0])
    top_middle = figure.add_subplot(grid[0, 1])
    top_text = figure.add_subplot(grid[0, 2:])
    bottom_axes = [figure.add_subplot(grid[1, index]) for index in range(4)]
    for axis, column, title in [(top_left, "nonsynonymous_mutation_burden", "Mutation burden"), (top_middle, "non_diploid_CNV_fraction", "Non-diploid CNV fraction")]:
        rho, p_value = spearmanr(genomic[column], genomic["VM_RNA_log2TPM"])
        axis.scatter(genomic[column], genomic["VM_RNA_log2TPM"], s=12, alpha=0.65, color=COLORS["grey"], edgecolor="none")
        axis.set_xlabel(title)
        axis.set_ylabel("RNA VM-like score")
        axis.text(0.02, 1.015, f"n={len(genomic)}\nρ={rho:.3f}\nP={p_value:.2e}", transform=axis.transAxes, va="bottom", fontsize=6.2,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 1.0, "pad": 0.25}, clip_on=False)
        axis.grid(color="#E9ECEF", linewidth=0.7)
    top_text.axis("off")
    top_text.text(0, 0.82, "Predefined genomic context", fontsize=9, fontweight="bold")
    top_text.text(0, 0.53, "KRAS mutation: higher score\nBAP1 mutation: lower score\nMYC/MCL1 gain: lower score", fontsize=8, linespacing=1.6)
    top_text.text(0, 0.12, "Association-only\nGenome-wide discrete-CNV scan: no BH-FDR-positive hit", fontsize=7, color="#8B0000", linespacing=1.6)
    event_panels = [("KRAS_mutant", "KRAS mutation", "FDR=0.0013"), ("BAP1_mutant", "BAP1 mutation", "FDR=0.0443"), ("MYC_gain", "MYC gain", "FDR=0.0029"), ("MCL1_gain", "MCL1 gain", "FDR=0.0029")]
    for axis, (column, title, fdr_label) in zip(bottom_axes, event_panels):
        wildtype = genomic.loc[~genomic[column], "VM_RNA_log2TPM"]
        altered = genomic.loc[genomic[column], "VM_RNA_log2TPM"]
        axis.boxplot([wildtype, altered], tick_labels=["WT/diploid", "Altered"], patch_artist=True, boxprops={"facecolor": "#D9EAF7"}, medianprops={"color": COLORS["dark"]})
        jitter_one = np.random.default_rng(17).normal(1, 0.035, len(wildtype))
        jitter_two = np.random.default_rng(29).normal(2, 0.035, len(altered))
        axis.scatter(jitter_one, wildtype, s=6, color=COLORS["grey"], alpha=0.45)
        axis.scatter(jitter_two, altered, s=6, color=COLORS["vm"], alpha=0.55)
        axis.set_title(title, fontsize=7.5, pad=30)
        axis.set_ylabel("RNA VM-like score")
        axis.text(0.03, 1.015, f"n altered={len(altered)}\n{fdr_label}", transform=axis.transAxes, va="bottom", fontsize=5.8,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 1.0, "pad": 0.25}, clip_on=False)
    panel_label(top_left, "A", x=-0.25, y=1.15)
    panel_label(top_middle, "B", x=-0.25, y=1.15)
    figure.text(0.045, 0.414, "C", fontsize=9, fontweight="bold", va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.15})
    figure.suptitle("Fig. 6 | Association-only genomic context of the VM-like epithelial state", x=0.01, y=1.01, ha="left", fontsize=10, fontweight="bold")
    figure.subplots_adjust(top=0.86, bottom=0.10, left=0.08, right=0.98)
    save_figure(figure, "Fig6_association_only_genomic_context")


def main():
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    setup_style()
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    figure_5()
    figure_6()


if __name__ == "__main__":
    main()
