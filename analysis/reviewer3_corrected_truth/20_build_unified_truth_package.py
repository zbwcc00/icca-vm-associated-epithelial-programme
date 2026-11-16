from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import RGBColor
from PIL import Image

ROOT = Path(r"D:\第六篇大论文")
R3 = ROOT / "revision_round1_reviewer3_20261116"
OLD = ROOT / "submission_package" / "discover_oncology_analysis_20261116_revised_round1_final_all_reviewers"
OUT = ROOT / "submission_package" / "discover_oncology_analysis_20261116_unified_truth_final"
WORK = R3 / "unification_work"
FIG = R3 / "figures_unified_truth"
RED = RGBColor(0xC0, 0x00, 0x00)
VM = "#CC79A7"
BLUE = "#0072B2"
TEAL = "#009E73"
ORANGE = "#E69F00"
GREY = "#B7C3CE"
DARK = "#263238"


def set_p(p, text, marked=False):
    for r in p.runs:
        r.text = ""
    run = p.runs[0] if p.runs else p.add_run()
    run.text = text
    if marked:
        run.font.color.rgb = RED


def add_after(p, text, marked=False, style=None):
    q = p._parent.add_paragraph(style=style)
    p._p.addnext(q._p)
    r = q.add_run(text)
    if marked:
        r.font.color.rgb = RED
    return q


def replace_docx_text(src, dst, replacements, marked=False):
    d = Document(src)
    for i, text in replacements.items():
        if i < len(d.paragraphs):
            set_p(d.paragraphs[i], text, marked)
    # Replace old terminology throughout body, leaving reviewer quotations in response files alone.
    for p in d.paragraphs:
        t = p.text
        if not t:
            continue
        u = t
        for a, b in {
            "final 30 groups": "corrected 26-group candidate set",
            "final 30": "corrected 26-group",
            "complete 28-group": "complete corrected 26-group",
            "complete 28 group": "complete corrected 26-group",
            "prespecified": "predefined",
            "Prespecified": "Predefined",
            "RNA–protein concordance": "RNA–protein correlation",
            "cross-layer concordance": "cross-layer association",
            "statistical concordance": "measurement agreement",
            "VM-associated": "VM-related",
            "VM-like": "VM-related",
        }.items():
            u = u.replace(a, b)
        if u != t:
            set_p(p, u, marked)
    d.save(dst)


def make_fig4():
    x = pd.read_csv(R3 / "analysis_10_corrected_bulk_candidate_audit/Table_U14_corrected_26_group_bulk_composition_audit.tsv", sep="\t")
    w = pd.read_csv(R3 / "analysis_10_corrected_bulk_candidate_audit/Table_U15_corrected_26_group_three_bulk_direction.tsv", sep="\t")
    keep = w.loc[w.positive_all_three].copy()
    order = keep.sort_values("GSE32225").signature_group.tolist()
    y = keep.set_index("signature_group").loc[order]
    names = y.external_source_RNA_members.str.replace(";", "; ", regex=False).tolist()
    fig = plt.figure(figsize=(11, 8.5), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=.55, wspace=.42)
    ax = fig.add_subplot(gs[0, :])
    ypos = np.arange(len(order))
    for j, (cohort, color) in enumerate([("GSE32225", BLUE), ("GSE107943", ORANGE), ("GSE26566", TEAL)]):
        ax.scatter(y[cohort], ypos + (j-1)*.18, s=28, c=color, label=cohort, zorder=3)
    ax.axvline(0, color="#808A93", lw=.8)
    ax.set_yticks(ypos); ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Composition-adjusted association coefficient (β)")
    ax.set_title("12 of 24 fully evaluable corrected groups were positive in all three cohorts", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=7); ax.grid(axis="x", color="#E9ECEF", lw=.6)
    ax.text(-.07, 1.08, "A", transform=ax.transAxes, fontweight="bold", fontsize=11)
    # Heatmap of all 26 group directions, keeping missingness explicit.
    allg = w.sort_values("signature_group")
    mat = allg[["GSE32225", "GSE107943", "GSE26566"]].to_numpy(float)
    h = fig.add_subplot(gs[1, 0])
    im = h.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    h.set_xticks(range(3)); h.set_xticklabels(["GSE32225", "GSE107943", "GSE26566"], rotation=35, ha="right", fontsize=6)
    h.set_yticks(range(len(allg))); h.set_yticklabels(allg.external_source_RNA_members.str.replace(";", "; ", regex=False), fontsize=5.3)
    h.set_title("Direction across the corrected 26-group family", fontsize=8, loc="left", fontweight="bold")
    fig.colorbar(im, ax=h, fraction=.045, pad=.03, label="β")
    h.text(-.18, 1.04, "B", transform=h.transAxes, fontweight="bold", fontsize=11)
    s = x.pivot(index="signature_group", columns="cohort", values="four_proxy_beta").reindex(allg.signature_group)
    q = fig.add_subplot(gs[1, 1])
    q.imshow(s.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    q.set_xticks(range(3)); q.set_xticklabels(["32225", "107943", "26566"], fontsize=6)
    q.set_yticks(range(len(s))); q.set_yticklabels(allg.external_source_RNA_members.str.replace(";", "; ", regex=False), fontsize=5.3)
    q.set_title("Sensitivity model coefficients", fontsize=8, loc="left", fontweight="bold")
    q.text(-.18, 1.04, "C", transform=q.transAxes, fontweight="bold", fontsize=11)
    fig.suptitle("Fig. 4 | Corrected 26-group bulk association audit", x=.01, ha="left", fontsize=11, fontweight="bold")
    fig.text(.5, .01, "Bulk covariance after expression-derived composition adjustment; this is directional association and not cellular-source or mechanistic validation.", ha="center", fontsize=6.7)
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "Fig4_corrected_26_group_bulk_audit.png", bbox_inches="tight", dpi=600)
    fig.savefig(FIG / "Fig4_corrected_26_group_bulk_audit.pdf", bbox_inches="tight")
    (OUT / "source_data_and_supplementary_tables").mkdir(parents=True, exist_ok=True)
    y.reset_index().to_csv(OUT / "source_data_and_supplementary_tables/Fig4_corrected_three_cohort_positive_groups.tsv", sep="\t", index=False)
    allg.to_csv(OUT / "source_data_and_supplementary_tables/Fig4_corrected_26_group_direction.tsv", sep="\t", index=False)
    x.to_csv(OUT / "source_data_and_supplementary_tables/Fig4_corrected_26_group_composition_audit.tsv", sep="\t", index=False)


def make_s15_s17_s18():
    FIG.mkdir(exist_ok=True)
    # S15: redundancy, anchor sensitivity, random calibration.
    red = pd.read_csv(R3 / "analysis_09_corrected_secact_redundancy/Table_U13_corrected_SecAct_signature_redundancy.tsv", sep="\t")
    rnd = pd.read_csv(R3 / "analysis_08_corrected_random_benchmark/Table_U12_corrected_detection_expression_random_benchmark.tsv", sep="\t")
    endo = pd.read_csv(R3 / "analysis_06_unified_corrected_sensitivities/Table_U10_corrected_independent_endothelial_exclusion.tsv", sep="\t")
    endo = endo.loc[endo.independent_endothelial_exclusion_quantile.eq(0.9)].copy()
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.7), dpi=300)
    ax[0].bar(red.set.astype(str), red.shannon_effective_rank, color=[GREY, BLUE, TEAL, VM, ORANGE])
    ax[0].set_ylabel("Shannon effective rank"); ax[0].tick_params(axis="x", rotation=35, labelsize=6)
    ax[0].set_title("Correlated hypothesis structure", loc="left", fontweight="bold", fontsize=8)
    ax[1].bar(rnd.cohort, rnd.empirical_two_sided_p, color=[BLUE, TEAL]); ax[1].axhline(.05, ls="--", color="#808A93")
    ax[1].set_ylabel("Empirical two-sided P"); ax[1].set_title("Detection/expression-matched null", loc="left", fontweight="bold", fontsize=8)
    ax[2].bar(endo.cohort, endo.fraction_original_high_score_retained, color=[BLUE, TEAL]); ax[2].set_ylim(0,1.05); ax[2].set_ylabel("Retained fraction")
    ax[2].set_title("Independent endothelial exclusion", loc="left", fontweight="bold", fontsize=8)
    for i,a in enumerate(ax): a.text(-.12,1.06,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=10)
    fig.suptitle("Fig. S15 | Corrected internal calibration and exclusion sensitivity", x=.01, ha="left", fontsize=10, fontweight="bold")
    fig.text(.5,.01,"Random benchmark is internal calibration because labels use the same six anchors; no panel is independent biological validation.",ha="center",fontsize=6)
    fig.savefig(FIG/"FigS15_corrected_internal_calibration.png",bbox_inches="tight",dpi=600); fig.savefig(FIG/"FigS15_corrected_internal_calibration.pdf",bbox_inches="tight"); plt.close(fig)
    # S17: corrected QC and scoring correlations.
    qc = pd.read_csv(R3 / "analysis_06_unified_corrected_sensitivities/Table_U5_corrected_patient_QC_summaries.tsv", sep="\t")
    cor = pd.read_csv(R3 / "analysis_06_unified_corrected_sensitivities/Table_U4_corrected_score_correlations.tsv", sep="\t")
    mit = pd.read_csv(R3 / "analysis_06_unified_corrected_sensitivities/Table_U7_corrected_mitochondrial_ceiling_sensitivity.tsv", sep="\t")
    fig, ax = plt.subplots(1,3,figsize=(11,3.7),dpi=300)
    qwide = qc.pivot_table(index=["cohort","patient_id"], columns=["score_group","metric"], values="median")
    for cohort,g in qwide.reset_index().groupby("cohort"):
        xcol=("high","n_counts"); ycol=("high","pct_mito")
        if xcol in g.columns and ycol in g.columns: ax[0].scatter(g[xcol], g[ycol], label=cohort)
    ax[0].set_xlabel("Median high-score UMI count"); ax[0].set_ylabel("Median high-score mitochondrial fraction"); ax[0].legend(frameon=False,fontsize=6)
    cor["comparison"] = cor["score_1"].str.replace("_full6", "", regex=False) + " vs " + cor["score_2"].str.replace("_full6", "", regex=False)
    piv=cor.pivot(index="cohort",columns="comparison",values="spearman_rho"); ax[1].imshow(piv.to_numpy(float),vmin=0,vmax=1,cmap="viridis",aspect="auto"); ax[1].set_xticks(range(len(piv.columns))); ax[1].set_xticklabels(piv.columns,rotation=35,ha="right",fontsize=6); ax[1].set_yticks(range(len(piv))); ax[1].set_yticklabels(piv.index,fontsize=6); ax[1].set_title("Score implementation agreement",loc="left",fontsize=8,fontweight="bold")
    mit2=mit.loc[mit.score_group.eq("all_malignant")].copy()
    for cohort,g in mit2.groupby("cohort"): ax[2].plot(g.mitochondrial_ceiling_percent,g.retained_fraction,marker="o",label=cohort)
    ax[2].set_xlabel("Mitochondrial ceiling"); ax[2].set_ylabel("Retained fraction"); ax[2].legend(frameon=False,fontsize=6)
    for i,a in enumerate(ax): a.text(-.12,1.06,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=10)
    fig.suptitle("Fig. S17 | Corrected QC, scoring and mitochondrial sensitivity",x=.01,ha="left",fontsize=10,fontweight="bold")
    fig.savefig(FIG/"FigS17_corrected_qc_score_sensitivity.png",bbox_inches="tight",dpi=600); fig.savefig(FIG/"FigS17_corrected_qc_score_sensitivity.pdf",bbox_inches="tight"); plt.close(fig)
    # S18: corrected comparator, random, SecAct adjustment.
    comp=pd.read_csv(R3/"analysis_06_unified_corrected_sensitivities/Table_U8_corrected_comparator_programmes.tsv",sep="\t")
    adj=pd.read_csv(R3/"analysis_07_corrected_q80_core_program_adjustment/Table_U11_corrected_q80_core_program_adjustment.tsv",sep="\t")
    fig,ax=plt.subplots(1,3,figsize=(11,3.7),dpi=300)
    for cohort,g in comp.groupby("cohort"): ax[0].bar(np.arange(5)+(0 if cohort=="GSE138709" else .35),g.hedges_g_high_vs_remainder,width=.35,label=cohort)
    ax[0].set_xticks(np.arange(5)+.175); ax[0].set_xticklabels(["EMT","Hypoxia","ECM","Migration","Endothelial"],rotation=35,ha="right",fontsize=6); ax[0].set_ylabel("Hedges g"); ax[0].legend(frameon=False,fontsize=6); ax[0].set_title("Anchor-excluded comparator programmes",loc="left",fontsize=8,fontweight="bold")
    ax[1].bar(rnd.cohort,rnd.six_anchor_paired_Hedges_g,color=[BLUE,TEAL]); ax[1].set_ylabel("Observed Hedges g"); ax[1].set_title("Matched random calibration",loc="left",fontsize=8,fontweight="bold")
    q=adj.query("candidate_set == 'corrected_q80_external_26'"); labels=q.adjustment.str.replace("_anchor_excluded", "", regex=False).str.replace("Hallmark_", "", regex=False).tolist();
    for cohort,g in q.groupby("cohort"): ax[2].bar(np.arange(len(g))+(0 if cohort=="GSE138709" else .4),g.spearman_zscore_unadjusted_vs_adjusted,width=.4,label=cohort)
    ax[2].set_xticks(np.arange(len(labels))+.2); ax[2].set_xticklabels(labels,rotation=35,ha="right",fontsize=6); ax[2].set_ylim(0,1.05); ax[2].set_ylabel("Unadjusted/adjusted ρ"); ax[2].legend(frameon=False,fontsize=6); ax[2].set_title("One-at-a-time SecAct overlap audit",loc="left",fontsize=8,fontweight="bold")
    for i,a in enumerate(ax): a.text(-.12,1.06,chr(65+i),transform=a.transAxes,fontweight="bold",fontsize=10)
    fig.suptitle("Fig. S18 | Corrected programme overlap, random calibration and SecAct adjustment",x=.01,ha="left",fontsize=10,fontweight="bold")
    fig.text(.5,.01,"These analyses quantify overlap and sensitivity; they do not establish VM specificity or independent validation.",ha="center",fontsize=6)
    fig.savefig(FIG/"FigS18_corrected_overlap_random_secact.png",bbox_inches="tight",dpi=600); fig.savefig(FIG/"FigS18_corrected_overlap_random_secact.pdf",bbox_inches="tight"); plt.close(fig)


def sync_tables_and_figures():
    sd = OUT / "source_data_and_supplementary_tables"
    s7 = OUT / "supplementary_tables_S7"
    # Replace label-dependent old tables with corrected truth tables.
    mapping = {
        "Supplementary_Table_S5f_SecAct_signature_redundancy.tsv": R3/"analysis_09_corrected_secact_redundancy/Table_U13_corrected_SecAct_signature_redundancy.tsv",
        "Supplementary_Table_S5g_SecAct_covariate_stability.tsv": R3/"analysis_07_corrected_q80_core_program_adjustment/Table_U11_corrected_q80_core_program_adjustment.tsv",
        "Supplementary_Table_S5j_bulk_composition_sensitivity.tsv": R3/"analysis_10_corrected_bulk_candidate_audit/Table_U14_corrected_26_group_bulk_composition_audit.tsv",
        "Supplementary_Table_S5k_independent_endothelial_exclusion.tsv": R3/"analysis_06_unified_corrected_sensitivities/Table_U10_corrected_independent_endothelial_exclusion.tsv",
        "Supplementary_Table_S6b1_patient_QC_summaries.tsv": R3/"analysis_06_unified_corrected_sensitivities/Table_U5_corrected_patient_QC_summaries.tsv",
        "Supplementary_Table_S6b2_patient_QC_tests.tsv": R3/"analysis_06_unified_corrected_sensitivities/Table_U6_corrected_patient_QC_tests.tsv",
        "Supplementary_Table_S6b3_mitochondrial_sensitivity.tsv": R3/"analysis_06_unified_corrected_sensitivities/Table_U7_corrected_mitochondrial_ceiling_sensitivity.tsv",
        "Supplementary_Table_S6c_score_correlations.tsv": R3/"analysis_06_unified_corrected_sensitivities/Table_U4_corrected_score_correlations.tsv",
        "Supplementary_Table_S6f_detection_expression_matched_random_benchmark.tsv": R3/"analysis_08_corrected_random_benchmark/Table_U12_corrected_detection_expression_random_benchmark.tsv",
    }
    for name, src in mapping.items():
        if src.exists(): shutil.copy2(src, sd/name)
    shutil.copy2(R3/"analysis_10_corrected_bulk_candidate_audit/Table_U14_corrected_26_group_bulk_composition_audit.tsv", s7/"Supplementary_Table_S7k_corrected_26_group_bulk_composition_audit.tsv")
    shutil.copy2(R3/"analysis_10_corrected_bulk_candidate_audit/Table_U15_corrected_26_group_three_bulk_direction.tsv", s7/"Supplementary_Table_S7l_corrected_26_group_three_bulk_direction.tsv")
    # Supplementary figures.
    supp = OUT / "supplementary_figures_S1-S19"
    for oldname, newbase in [("FigS15_robustness_calibration_composite", "FigS15_corrected_internal_calibration"), ("FigS17_cellular_origin_and_score_audits", "FigS17_corrected_qc_score_sensitivity"), ("FigS18_program_random_F3_protein_CNV_audits", "FigS18_corrected_overlap_random_secact")]:
        for ext in ("pdf", "png"):
            src=FIG/(newbase+"."+ext)
            if src.exists(): shutil.copy2(src, supp/(oldname+"."+ext))
    # Main figure 4 is embedded in manuscript; keep standalone corrected assets too.
    for ext in ("pdf", "png"):
        shutil.copy2(FIG/("Fig4_corrected_26_group_bulk_audit."+ext), OUT/("Fig4_corrected_26_group_bulk_audit."+ext))


def build_docs():
    clean_repl = {
        9: "Results: After complete six-gene reconstruction and endothelial-control filtering, q80 operational high-score cells comprised 1,331/12,447 cells in GSE138709 (10.69%; 4/4 patients), 220/2,840 in GSE181878 (7.75%; one-patient feasibility) and 277/1,888 in GSE125449 (14.67%; 9/10 patients). At q90, GSE138709 retained 633 cells (5.09%; 4/4 patients) and GSE125449 retained 122 cells (6.46%; 7/10 patients). De novo leave-one-anchor reclassification was generally stable for CDH5 and MMP9, moderately sensitive to EPHA2 and MMP2, and most sensitive to LAMC2. The corrected q80 SecAct funnel was 657→268→58→26; the q90 analysis retained 244, 56 and 26 groups, with receiver z-score correlations of 0.982 and 0.967 between thresholds. In 208 matched tumours, the RNA and total-protein programme scores were correlated (Spearman rho = 0.709). Phosphopeptide and genomic associations were reported with sample sizes, confidence intervals, nominal P values and within-family BH values.",
        52: "For each cohort, operational malignant epithelial cells were aggregated by patient or library and, where applicable, high-score versus remainder. Raw UMI counts were summarised as pseudobulk profiles and analysed with edgeR [22]. At q80, GSE138709 contributed four eligible paired libraries from three patients (eight pseudobulks; design ~ sample_id + score_group; minimum high-score count, five cells), and GSE125449 contributed five paired libraries from five patients (ten pseudobulks; design ~ library_id + score_group). At q90, the same rule yielded four GSE138709 libraries from three patients and four GSE125449 libraries from four patients. These paired contrasts provide within-library programme associations and do not estimate independent patient-level effects.",
        116: "We re-ran the corrected 26-group family in tumour-only bulk samples from GSE32225 (n = 149), GSE107943 (n = 30) and GSE26566 (n = 104). Of 24 groups evaluable in all three cohorts, 12 had positive composition-adjusted coefficients in all three cohorts and 17 had the same direction across the three cohorts. KCP and PDGFA were not fully evaluable across platforms because member expression was unavailable in one or more cohorts. The complete corrected family yielded 25 evaluable groups in GSE32225, 26 in GSE107943 and 24 in GSE26566; positive adjusted 95% confidence intervals were observed for 11, 3 and 6 groups, respectively. MMP2-removal estimates were recomputed in parallel. These models reduce measured composition confounding but cannot account fully for tumour purity, malignant-cell abundance, subtype, platform or shared-matrix covariance. The results therefore show directional bulk association only, not the same cellular origin or mechanism across cohorts (Fig. 4 and corrected bulk source tables).",
        133: "The treatment-response audit included nine responders and nine non-responders in GSE255058. Neither the VM-related score nor PLAU/SERPINE1, IL1A or S100A8/S100A9 differed after exact permutation, bootstrap or myeloid-adjusted HC3 analyses. In GSE244807, an unadjusted OS association in 230 evaluable samples (HR, 1.247 per standard deviation; P = 0.0109) was attenuated after tissue-type adjustment (HR, 0.979; P = 0.835) and joint adjustment (HR, 0.958; P = 0.676). Tissue-stratified estimates were also non-robust. In the BTC spatial audit (six patients and 12 sections), patient-level VM-related and myeloid-state score correlations were small and unsupported by corrected sign tests. Finally, cross-cohort ExtraTrees counterfactuals and scTenifoldKnk perturbations did not nominate a dual-algorithm upstream target. CRISPR analyses yielded no FDR-positive dependency. The sole PRISM signal, SF-11, remained a single-dose in-vitro association without clinical interpretability.",
        159: "SecAct is used here only to prioritise hypotheses. The corrected 26 groups show substantial target-signature redundancy, with effective rank 4.08 and 34.46% of pairwise correlations at |r| ≥ 0.7. One-at-a-time adjustment of the corrected q80 family gave z-score correlations of 0.448, 0.975 and 0.813 for EMT, hypoxia and ECM in GSE138709, and 0.804, 0.738 and 0.731 in GSE125449. Receiver-side activity, source-compartment RNA and bulk covariance are separate observations; none measures secretion, receptor engagement, spatial proximity or functional communication. A defensible validation sequence is to quantify compartment-resolved protein, test spatial proximity and receptor activity, and then perturb candidate pathways in an organoid–macrophage system.",
        208: "Ethical approval",
        209: "Not applicable. This study analysed only publicly available, de-identified datasets and involved no recruitment, intervention, prospective sample collection or access to identifiable information. The study followed the ethical principles of the Declaration of Helsinki as applicable to secondary analyses of publicly available human data and the governance requirements of the originating studies. No separate approval committee or institutional review board approval was required for this secondary analysis.",
        210: "Consent to participate",
        211: "Not applicable. No participant recruitment or direct interaction occurred in this secondary analysis of de-identified public data.",
        212: "Consent to publish",
        213: "Not applicable. No identifiable individual-level data, photographs or case details are presented.",
    }
    # Existing paragraphs after 213 are shifted only if the declaration insertion is not performed.
    for marked in (False, True):
        src = OLD / ("02_Manuscript_Revised_Red_Marked.docx" if marked else "01_Manuscript_Revised_Clean.docx")
        dst = OUT / ("02_Manuscript_Revised_Red_Marked.docx" if marked else "01_Manuscript_Revised_Clean.docx")
        replace_docx_text(src, dst, clean_repl, marked)
    # Fix declaration structure after global rewrite by replacing any residual declaration wording.
    for fn, marked in [("01_Manuscript_Revised_Clean.docx",False),("02_Manuscript_Revised_Red_Marked.docx",True)]:
        d=Document(OUT/fn)
        for p in d.paragraphs:
            if p.text == "Consent for publication": set_p(p,"Consent to publish",marked)
            if "Ethics approval and consent to participate" in p.text: set_p(p,"Ethical approval",marked)
        d.save(OUT/fn)


def update_response_doc(src, dst, replacement_map, title=None):
    d=Document(src)
    if title:
        for p in d.paragraphs:
            if p.text.startswith("Manuscript title:"):
                set_p(p, "Manuscript title: "+title, False)
    for i,t in replacement_map.items():
        if i < len(d.paragraphs): set_p(d.paragraphs[i],t,False)
    for p in d.paragraphs:
        t=p.text; u=t
        for a,b in {"prespecified":"predefined","Prespecified":"Predefined","final-30":"corrected 26-group","final 30":"corrected 26-group","657→249→66→30":"657→268→58→26 (q80) and 657→244→56→26 (q90)","P = 0.029":"P = 0.011","P = 0.222":"P = 0.155","concordance":"correlation","vasculogenic-mimicry-associated":"VM-related","VM-associated":"VM-related","VM-like":"VM-related","patient-covered state":"patient-covered continuous programme"}.items(): u=u.replace(a,b)
        if u!=t: set_p(p,u,False)
    d.save(dst)


def build_responses():
    title="Cross-cohort definition and multi-omic characterization of a six-gene VM-related epithelial transcriptional programme in intrahepatic cholangiocarcinoma"
    # Reviewer 1 response paragraphs are 6, 12, 18, 24, 30, 36, 42, 48, 50 in the existing file.
    r1={
      2:"We thank Reviewer 1 for identifying the circularity and patient-level limitations. We rebuilt all label-dependent analyses from the complete six-gene raw-UMI reconstruction and replaced the old shortlist-derived results. The corrected 26-group bulk audit and all corrected source tables are included. Claims remain bounded to a continuous VM-related programme and association-level calibration.",
      6:"We agree that the random benchmark is not independent validation. The corrected detection- and expression-matched benchmark is now reported only as internal calibration, with empirical P = 0.011 in GSE138709 and P = 0.155 in GSE125449. Anchor-excluded comparator programmes and de novo leave-one-anchor analyses quantify overlap and dependence, but no morphology or perfusion endpoint was available.",
      12:"We agree that the anchors are not VM-specific markers. Corrected anchor-excluded Hedges g values in GSE138709 were 1.163 for EMT, 0.466 for hypoxia, 1.188 for ECM, 0.417 for migration and 0.686 for endothelial identity; corresponding GSE125449 values were 0.236, 0.162, 0.391, 0.213 and 0.488. We therefore use continuous VM-related programme terminology and do not claim molecular specificity.",
      18:"We agree that dataset recurrence is not patient-level replication. The corrected leave-one-patient-out Jaccard ranges were 0.726–0.942 in GSE138709 and 0.628–1.000 in GSE125449, with corrected pseudobulk eligibility of 4 libraries/3 patients at q80 and q90 in GSE138709, 5/5 at q80 and 4/4 at q90 in GSE125449. GSE181878 remains single-patient feasibility evidence.",
      24:"We corrected the role of the endothelial-control filter and recomputed the independent exclusion using the complete six-gene scores. At the q90 score cutoff, the independent endothelial-identity exclusion retained 1,038/1,331 GSE138709 high-score cells (77.99%, all four patients) and 233/277 GSE125449 cells (84.12%, nine patients). CNV-reference sensitivities remain separate operational classification audits. These results reduce, but cannot eliminate, ambient or admixture explanations.",
      30:"We recalculated the threshold grid after restoring all six anchors. q80/q90 prevalence, patient coverage and eligibility are reported in the corrected tables. SecAct receiver z-score correlations between q80 and q90 were 0.982 in GSE138709 and 0.967 in GSE125449, with 233 nominal groups overlapping (Jaccard 0.835). The three scoring implementations are correlated summaries of one anchor set, not independent evidence.",
      36:"We replaced the old 249/66/30 and final-30 results. The corrected funnel is 657→268→58→26 at q80 and 657→244→56→26 at q90. The corrected 26-group family has effective rank 4.08, with 34.46% of pairwise correlations at |r| ≥ 0.7. One-at-a-time adjustment correlations for EMT, hypoxia and ECM were 0.448, 0.975 and 0.813 in GSE138709 and 0.804, 0.738 and 0.731 in GSE125449. These are overlap audits, not evidence of independent SecAct mechanisms.",
      42:"We retain the source-RNA result only as non-exclusive compatibility. PLAU is detectable in myeloid and stromal compartments, while SERPINE1 is also prominent in stromal and malignant epithelial compartments. GSE201425 has two patients and no epithelial high-score/remainder contrast. We therefore do not assign a macrophage-specific signalling axis.",
      48:"We re-ran the complete corrected 26-group family in all three bulk cohorts. Of 24 groups evaluable in all three, 12 were positive in all three and 17 had the same direction. The cohort-specific complete-family counts were 25, 26 and 24 evaluable groups, with positive adjusted 95% confidence intervals for 11, 3 and 6 groups in GSE32225, GSE107943 and GSE26566, respectively. KCP and PDGFA were not fully evaluable across platforms. MMP2-removal estimates were also recomputed. These remain directional bulk associations, not cellular-source or mechanistic validation.",
      50:"We tested all corrected 26 candidate groups in GSE32225 (n = 149), GSE107943 (n = 30) and GSE26566 (n = 104). Twelve of 24 fully evaluable groups were positive in all three cohorts and 17 showed the same direction. Composition-proxy and MMP2-removal estimates are supplied in the corrected bulk tables. Histological purity, malignant-cell abundance and harmonised subtype covariates remain unavailable across the three cohorts, so this point is addressed with a bounded directional audit rather than claimed as independent validation.",
    }
    update_response_doc(OLD/"03_Response_to_Reviewer_1.docx", OUT/"03_Response_to_Reviewer_1.docx", r1)
    r2={
      1:"Manuscript title: "+title,
      5:"We agree that the endothelial-control rule cannot quantify diffuse ambient RNA. We therefore retained lineage-resolved detection and expression, corrected patient-level QC and mitochondrial-ceiling audits, independent endothelial exclusion, and DepMap biliary-line expression. Formal SoupX, CellBender or DecontX correction was not estimable because the retained public objects contain filtered cells without empty-droplet or unfiltered background matrices. The corrected analyses assess compatibility and reduce obvious doublet or composition explanations, but do not prove tumour-cell-exclusive origin.",
      8:"We agree that the dataset-specific mitochondrial ceilings are permissive. We now call the 15% and 20% analyses predefined sensitivity analyses. In the corrected cohorts, the high-score fractions retained under these ceilings are reported in the corrected Table U7, with GSE138709 retaining 79.3% and 90.8% of q80 high-score cells and GSE125449 retaining 86.6% and 91.0%, respectively. These analyses do not prove that stressed cells were removed.",
      11:"We agree. UCell, AUCell and the module score are correlated implementations of the same six-gene score. Corrected pairwise implementation correlations are reported in Table U4 and Fig. S17C. Their agreement is an operational robustness check, not independent biological validation.",
      14:"F3 was observed during transcript discovery and was later excluded from soluble-ligand interpretation because it is absent from the SecAct v1.1 catalogue. The corrected catalogue-defined funnels are 657→268→58→26 at q80 and 657→244→56→26 at q90. F3 inclusion therefore cannot alter the catalogue-defined funnel; it remains a transcript-level state-support audit.",
      17:"We replaced all old final-30 adjustment values with corrected 26-group one-at-a-time audits. For EMT, hypoxia and ECM, q80 unadjusted-versus-adjusted correlations were 0.448, 0.975 and 0.813 in GSE138709 and 0.804, 0.738 and 0.731 in GSE125449. Coverage partial correlation remains inappropriate because coverage is threshold-derived from the same score and the patient denominators are small. RNA–protein is reported as a cross-layer correlation, not an independent endpoint.",
      20:"We repeated the detection- and expression-matched benchmark with 1,000 fixed-seed random six-gene sets. The corrected empirical two-sided P values were 0.011 for GSE138709 and 0.155 for GSE125449. Because the labels are defined by the same anchors, this remains internal calibration rather than independent validation.",
      23:"All six proteins were quantified in all 208 strictly matched OEP001105 tumours without imputation. The corrected RNA–protein programme correlation was rho = 0.709, P = 4.17 × 10−33, and the endothelial-composition partial correlation was rho = 0.688, P = 1.76 × 10−30. This is a within-project cross-layer association, not independent replication or formal purity adjustment.",
      29:"We harmonised the terminology as continuous six-gene VM-related epithelial programme, operational upper-tail label, external source-RNA compatibility audit, fixed-centroid cell-identity transfer, operational malignant epithelial cells and two-algorithm virtual perturbation. Definitions are supplied at first use.",
      41:"We audited all denominators after the raw-UMI correction. GSE138709 contributed four patients overall and 3 patients to paired q80/q90 pseudobulk eligibility; GSE125449 contributed 5 patients at q80 and 4 at q90; GSE201425 had 2 patients; OEP001105 had 208 matched RNA–protein tumours and 249 RNA–WES/CNV tumours with 204 complete genomic feature sets; the spatial audit had 6 BTC patients; and the treatment audit had 9 responders and 9 non-responders. Panel-specific source tables report the remaining sample-level denominators.",
    }
    update_response_doc(OLD/"04_Response_to_Reviewer_2.docx", OUT/"04_Response_to_Reviewer_2.docx", r2, title=title)
    # Reviewer 3 gets terminology and one denominator correction, while preserving its comments.
    update_response_doc(OLD/"05_Response_to_Reviewer_3.docx", OUT/"05_Response_to_Reviewer_3.docx", {2:"We thank Reviewer 3 for identifying the missing raw-UMI anchors. We now treat the complete six-gene reconstruction as the sole source of truth and have propagated it through scoring, labels, pseudobulk, SecAct, bulk audits, figures, source tables and the Reviewer 1/2 responses."})


def package():
    if OUT.exists(): shutil.rmtree(OUT)
    shutil.copytree(OLD, OUT)
    make_fig4(); make_s15_s17_s18(); build_docs(); build_responses(); sync_tables_and_figures()
    # Copy corrected figure source and analysis archive into package.
    arch=OUT/"reviewer3_analysis_archive"/"unified_truth_outputs"; arch.mkdir(parents=True,exist_ok=True)
    for p in [R3/"analysis_10_corrected_bulk_candidate_audit", R3/"analysis_06_unified_corrected_sensitivities", R3/"analysis_07_corrected_q80_core_program_adjustment", R3/"analysis_08_corrected_random_benchmark", R3/"analysis_09_corrected_secact_redundancy", FIG]:
        if p.exists(): shutil.copytree(p,arch/p.name)
    manifest={
      "source_of_truth":"Reviewer 3 complete six-gene raw-UMI reconstruction and corrected q80/q90 analyses",
      "truth_root":str(R3), "package_root":str(OUT), "corrected_secact_funnel":{"q80":[657,268,58,26],"q90":[657,244,56,26]},
      "corrected_bulk":{"fully_evaluable_groups":24,"positive_all_three":12,"same_direction_all_three":17,"cohort_evaluable":{"GSE32225":25,"GSE107943":26,"GSE26566":24}},
      "pseudobulk_eligibility":{"GSE138709":{"q80":"4 libraries/3 patients","q90":"4 libraries/3 patients"},"GSE125449":{"q80":"5 libraries/5 patients","q90":"4 libraries/4 patients"}},
      "superseded":{"old_shortlists":["657→249→66→30","final-30","complete 28-group audit"],"old_results":["old random P values","old final-30 adjustment values","old 28-group bulk counts","old endothelial-exclusion counts"]},
      "boundaries":["No morphology/perfusion claim","No secretion or communication claim","No independent clinical validation","Bulk analyses are directional covariance audits"]
    }
    (OUT/"10_unified_truth_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    (OUT/"00_README_Revision_Package_CN.md").write_text("# Unified truth revision package\n\n本包以 Reviewer 3 的完整六基因 raw-UMI 重算为唯一真值。旧 30 组、旧 28 组 bulk、旧 random P 值和旧标签派生表均已撤回或替换。\n\n- 主稿 clean 与 marked 已同步。\n- Reviewer 1/2 回复已按纠正后的结果重写。\n- Fig. 4、S15、S17、S18 及其 source tables 已重建。\n- corrected bulk tables U14/U15 与统一真值 manifest 随包提供。\n- 任何结果均不应解释为形态学 VM、分泌、通讯、因果、临床预测或独立验证。\n",encoding="utf-8")
    return OUT


if __name__ == "__main__":
    out=package(); print(out)
