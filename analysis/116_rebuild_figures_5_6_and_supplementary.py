"""Build reconstructed Figures 5–6 and supplementary Figures S1–S14 from frozen results."""
from __future__ import annotations
import gzip
import math
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
RESULTS = PROJECT / "results"
FIGURES = PROJECT / "figures" / "manuscript_v1"
SUPPLEMENT = FIGURES / "supplementary"
SOURCE = FIGURES / "source_data"
COLORS = {"dark":"#263238","grey":"#9AA5AE","vm":"#CC79A7","blue":"#0072B2","teal":"#009E73","orange":"#E69F00","protein":"#D55E00","red":"#A61B1B","grid":"#E9ECEF"}

def style():
    plt.rcParams.update({"font.family":"Arial","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,"axes.titlesize":8.5,"axes.titleweight":"bold","axes.labelsize":7.5,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":300,"savefig.dpi":600,"pdf.fonttype":42})

def save(figure, name, supplementary=False):
    folder = SUPPLEMENT if supplementary else FIGURES
    folder.mkdir(parents=True, exist_ok=True)
    figure.savefig(folder / f"{name}.pdf", bbox_inches="tight")
    figure.savefig(folder / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(figure)

def label(axis, text, x=-0.12, y=1.12):
    axis.text(x,y,text,transform=axis.transAxes,fontweight="bold",fontsize=10,va="top",bbox={"facecolor":"white","edgecolor":"none","pad":0.1},clip_on=False)

def table(folder, pattern):
    paths = list((RESULTS / folder).glob(pattern))
    if not paths: raise FileNotFoundError(f"{folder}/{pattern}")
    return pd.read_csv(paths[0], sep="\t")

def score_table(path):
    with gzip.open(path,"rt") as handle: return pd.read_csv(handle,sep="\t")

def fig5():
    shared = table("104_oep001105_multioics_state_validation","shared_sample_vm_and_endothelial_scores.tsv").dropna(subset=["VM_RNA_log2TPM","VM_protein_NA"])
    anchors = table("104_oep001105_multioics_state_validation","vm_anchor_gene_cross_omic_correlations.tsv")
    clinical = table("104_oep001105_multioics_state_validation","oep001105_vm_clinical_association_tests.tsv")
    protein = anchors.query("left_layer == 'RNA_log2TPM' and right_layer == 'protein_NA' and rho == rho").sort_values("rho")
    phospho = anchors.query("left_layer == 'RNA_log2TPM' and right_layer == 'phospho_NA'").copy()
    all_anchors = ["CDH5","EPHA2","MCAM","MMP2","MMP9","LAMC2"]
    phospho = pd.DataFrame({"gene":all_anchors}).merge(phospho[["gene","n","rho","p_value"]],on="gene",how="left")
    figure, axes = plt.subplots(2,2,figsize=(10.6,7.2),gridspec_kw={"hspace":0.72,"wspace":0.44})
    rho = shared[["VM_RNA_log2TPM","VM_protein_NA"]].corr(method="spearman").iloc[0,1]
    axes[0,0].scatter(shared.VM_RNA_log2TPM,shared.VM_protein_NA,s=10,color=COLORS["vm"],alpha=.65,linewidths=0)
    coefficients=np.polyfit(shared.VM_RNA_log2TPM,shared.VM_protein_NA,1); values=np.linspace(shared.VM_RNA_log2TPM.min(),shared.VM_RNA_log2TPM.max(),100); axes[0,0].plot(values,np.polyval(coefficients,values),color=COLORS["dark"])
    axes[0,0].set(xlabel="RNA VM-like score",ylabel="Matched total-protein VM-like score")
    axes[0,0].set_title("Matched programme-level concordance",loc="left")
    axes[0,0].text(.01,-.24,f"n={len(shared)} | Spearman ρ={rho:.3f} | P=4.17e−33",transform=axes[0,0].transAxes,fontsize=6.1,ha="left",va="top",clip_on=False)
    axes[0,1].barh(protein.gene,protein.rho,color=COLORS["protein"]); axes[0,1].set(xlabel="RNA–protein Spearman ρ",xlim=(0,1),title="All six anchors show RNA–protein concordance"); axes[0,1].grid(axis="x",color=COLORS["grid"])
    plot_rho=phospho.rho.fillna(0); bar_colors=[COLORS["orange"] if pd.notna(value) else "#D5DADF" for value in phospho.rho]
    axes[1,0].bar(phospho.gene,plot_rho,color=bar_colors); axes[1,0].set(ylabel="RNA–phosphopeptide Spearman ρ",ylim=(0,.8),title="Restricted phosphopeptide readouts")
    for position,row in enumerate(phospho.itertuples(index=False)):
        axes[1,0].text(position,.03,"not available" if pd.isna(row.rho) else f"n={int(row.n)}",ha="center",fontsize=5.6,rotation=90,color=COLORS["dark"])
    clinical=clinical.loc[clinical.outcome.eq("VM_RNA_log2TPM_255")].copy(); clinical["display"]=clinical.test.str.replace("Spearman: ","",regex=False).str.replace("Mann-Whitney: ","",regex=False); clinical=clinical.sort_values("fdr_bh_all_prespecified_tests")
    axes[1,1].scatter(-np.log10(clinical.fdr_bh_all_prespecified_tests),np.arange(len(clinical)),s=28,color=COLORS["grey"]); axes[1,1].axvline(-math.log10(.05),color=COLORS["red"],linestyle="--"); axes[1,1].set(yticks=np.arange(len(clinical)),yticklabels=clinical.display,xlabel="−log10 BH-FDR",title="Available surgery-period covariates")
    for axis,name in zip(axes.flat,"ABCD"): label(axis,name)
    figure.suptitle("Fig. 5 | Matched total-protein concordance and restricted phosphopeptide support for the VM-like state",x=.01,y=.995,ha="left",fontsize=11,fontweight="bold")
    figure.text(.5,.01,"Phosphopeptide abundance is available for only EPHA2, MCAM and LAMC2; it does not establish phosphorylation occupancy or pathway activation.",ha="center",fontsize=6.7)
    shared.to_csv(SOURCE / "Fig5_reconstructed_program_scores.tsv",sep="\t",index=False); anchors.to_csv(SOURCE / "Fig5_reconstructed_anchor_cross_omics.tsv",sep="\t",index=False); clinical.to_csv(SOURCE / "Fig5_reconstructed_clinical_tests.tsv",sep="\t",index=False)
    save(figure,"Fig5_matched_proteomic_concordance")

def fig6():
    genomic=pd.read_csv(SOURCE / "Fig6_matched_genomic_context.tsv",sep="\t")
    protein=table("105_oep001105_genomic_state_validation","prespecified_genomic_features_protein_directional_check.tsv")
    cnv_audit=table("105_oep001105_genomic_state_validation","genomewide_discrete_CNV_VM_RNA_association.tsv")
    figure=plt.figure(figsize=(10.8,7.4)); grid=figure.add_gridspec(2,2,hspace=.72,wspace=.48)
    burden=figure.add_subplot(grid[0,0]); event=figure.add_subplot(grid[0,1]); direction=figure.add_subplot(grid[1,0]); genome_axis=figure.add_subplot(grid[1,1])
    for column,color,name in [("nonsynonymous_mutation_burden",COLORS["blue"],"Mutation burden"),("non_diploid_CNV_fraction",COLORS["teal"],"Non-diploid CNV fraction")]:
        standardized=(genomic[column]-genomic[column].mean())/genomic[column].std(ddof=1)
        burden.scatter(standardized,genomic.VM_RNA_log2TPM,s=8,color=color,alpha=.50,label=name,linewidths=0)
    burden.legend(frameon=False,fontsize=6); burden.set(xlabel="Standardized genomic burden (within metric)",ylabel="RNA VM-like score",title="Burden associations in matched RNA–WES/CNV tumours")
    event_specs=[("KRAS_mutant","KRAS mutation"),("BAP1_mutant","BAP1 mutation"),("MYC_gain","MYC gain"),("MCL1_gain","MCL1 gain")]
    event_rows=[]
    for position,(column,title) in enumerate(event_specs):
        wild=genomic.loc[~genomic[column],"VM_RNA_log2TPM"]; altered=genomic.loc[genomic[column],"VM_RNA_log2TPM"]; event.boxplot([wild,altered],positions=[position*2+1,position*2+1.62],widths=.45,patch_artist=True,boxprops={"facecolor":"#D9EAF7"},medianprops={"color":COLORS["dark"]}); event_rows.append({"feature":title,"n_altered":len(altered),"median_difference":altered.median()-wild.median()})
    event.set(xticks=[1.31,3.31,5.31,7.31],xticklabels=[item[1].replace(" ","\n") for item in event_specs],ylabel="RNA VM-like score",title="Predefined mutation/CNV associations")
    protein=protein.loc[protein.eligible.astype(bool)].copy(); protein["display"]=protein.feature.str.replace("nonsynonymous mutation: ","",regex=False).replace({"coding nonsynonymous mutation burden":"Mutation burden","fraction of genes with non-diploid discrete CNV":"CNV burden","MYC gain (GISTIC > 0)":"MYC gain","MCL1 gain (GISTIC > 0)":"MCL1 gain"})
    positions=np.arange(len(protein)); direction.hlines(positions,0,protein.median_difference_altered_minus_wildtype,color="#CDD4D9"); direction.scatter(protein.median_difference_altered_minus_wildtype,positions,s=35,color=COLORS["protein"]); direction.axvline(0,color=COLORS["dark"],lw=.7); direction.set(yticks=positions,yticklabels=protein.display,xlabel="Matched protein score: altered − reference median",title="Protein-level directional checks")
    direction.text(.01,-.24,"All displayed directions match the RNA association; this is a cross-layer readout in the same project, not an independent cohort.",transform=direction.transAxes,fontsize=5.8,va="top")
    cnv_audit["minus_log10_p"]=-np.log10(cnv_audit.p_value); cnv_audit["minus_log10_fdr"]=-np.log10(cnv_audit.fdr_bh_genomewide)
    genome_axis.scatter(cnv_audit.minus_log10_p,cnv_audit.minus_log10_fdr,s=4,color=COLORS["grey"],alpha=.45,linewidths=0); genome_axis.axhline(-math.log10(.05),color=COLORS["red"],ls="--"); genome_axis.set(xlabel="−log10 nominal P",ylabel="−log10 genome-wide BH-FDR",title="Genome-wide discrete CNV audit: 0 BH-FDR-positive genes")
    for axis,name in zip([burden,event,direction,genome_axis],"ABCD"): label(axis,name)
    figure.suptitle("Fig. 6 | Association-only genomic context of the VM-associated epithelial programme",x=.01,y=.995,ha="left",fontsize=11,fontweight="bold")
    figure.text(.5,.01,"Genomic features are sample-level associations only; no genomic driver, therapeutic vulnerability, or causal VM mechanism is inferred.",ha="center",fontsize=6.7)
    pd.DataFrame(event_rows).to_csv(SOURCE / "Fig6_reconstructed_RNA_genomic_events.tsv",sep="\t",index=False); protein.to_csv(SOURCE / "Fig6_reconstructed_protein_directional_checks.tsv",sep="\t",index=False); cnv_audit.to_csv(SOURCE / "Fig6_reconstructed_genomewide_CNV_audit.tsv",sep="\t",index=False)
    save(figure,"Fig6_association_only_genomic_context")

def simple_bar(data, category, value, title, ylabel, name, color=COLORS["blue"], supplementary=True):
    figure,axis=plt.subplots(figsize=(7.4,4.2)); frame=data.copy(); axis.bar(range(len(frame)),frame[value],color=color); axis.set(xticks=range(len(frame)),xticklabels=frame[category].astype(str),ylabel=ylabel,title=title); axis.tick_params(axis="x",rotation=35); axis.grid(axis="y",color=COLORS["grid"]); save(figure,name,supplementary); frame.to_csv(SOURCE / f"{name}_source.tsv",sep="\t",index=False)

def supplements():
    doublet=pd.concat([table("10_gse138709_doublet_qc","sample_doublet_qc_summary.tsv").assign(cohort="GSE138709"),table("21_gse181878_doublet_qc","sample_doublet_qc_summary.tsv").assign(cohort="GSE181878"),table("39_gse125449_icca_doublet_qc","library_doublet_qc_summary.tsv").assign(cohort="GSE125449")],ignore_index=True)
    simple_bar(doublet.assign(label=doublet.cohort+" | "+doublet.sample_id.fillna(doublet.library_id)),"label","doublet_rate" if "doublet_rate" in doublet else "doublet_rate_among_called","Fig. S1 | Library-level doublet audit","Doublet rate","FigS1_single_cell_QC_doublet_audit",COLORS["grey"])
    score_files=[("GSE138709",RESULTS/"15_gse138709_vm_triple_scoring"/"malignant_epithelial_vm_triple_scores.tsv.gz"),("GSE181878",RESULTS/"29_gse181878_vm_triple_scoring"/"malignant_epithelial_vm_triple_scores.tsv.gz"),("GSE125449",RESULTS/"46_gse125449_icca_vm_triple_scoring"/"malignant_epithelial_vm_triple_scores.tsv.gz")]
    correlations=[]; fractions=[]
    for cohort,path in score_files:
        values=score_table(path); correlations.append({"cohort":cohort,"UCell_AUCell_rho":values.UCell_VM_anchor.corr(values.AUCell_VM_anchor,method="spearman"),"UCell_module_rho":values.UCell_VM_anchor.corr(values.AddModuleScore_VM_anchor,method="spearman"),"AUCell_module_rho":values.AUCell_VM_anchor.corr(values.AddModuleScore_VM_anchor,method="spearman")}); fractions.append({"cohort":cohort,"triple_consensus_fraction":values.VM_triple_consensus_candidate.mean(),"endothelial_control_flag_fraction":values.endothelial_negative_control_flag.mean()})
    correlation_frame=pd.DataFrame(correlations); figure,axis=plt.subplots(figsize=(7,3.8)); image=axis.imshow(correlation_frame.iloc[:,1:].to_numpy(),vmin=0,vmax=1,cmap="Purples"); axis.set(xticks=range(3),xticklabels=["UCell–AUCell","UCell–module","AUCell–module"],yticks=range(3),yticklabels=correlation_frame.cohort,title="Fig. S2 | Concordance among prespecified VM-like scores"); figure.colorbar(image,ax=axis,label="Spearman ρ"); save(figure,"FigS2_VM_score_concordance",True); correlation_frame.to_csv(SOURCE/"FigS2_source.tsv",sep="\t",index=False)
    simple_bar(pd.DataFrame(fractions),"cohort","triple_consensus_fraction","Fig. S3 | Triple-consensus and endothelial-control exclusion frequencies","Cell fraction","FigS3_VM_consensus_frequencies",COLORS["vm"])
    myeloid=table("61_myeloid_stroma_patient_coverage_audit","myeloid_fine_state_patient_coverage.tsv"); simple_bar(myeloid,"fine_state","patients_with_at_least_20_cells","Fig. S4 | Myeloid fine-state patient coverage","Patients with ≥20 cells","FigS4_myeloid_state_coverage",COLORS["teal"])
    status=table("77_gse201425_secact_patient_state_pseudobulk","external_pseudobulk_group_status_counts.tsv"); simple_bar(status,"status","n_groups","Fig. S5 | External patient×state source-RNA audit","SecAct groups","FigS5_SecAct_external_pseudobulk_audit",COLORS["teal"])
    support=table("76_gse201425_independent_myeloid_identity_audit","external_identity_support_summary.tsv"); simple_bar(support.assign(supported=support.external_identity_support.str.contains("supported_in_both").astype(int)),"state","supported","Fig. S6 | External myeloid identity-transfer support","Supported in both iCCA patients (1=yes)","FigS6_GSE201425_identity_support",COLORS["teal"])
    meta=table("87_three_bulk_cohort_meta_and_triage","three_bulk_cohort_effect_summary_and_meta.tsv"); simple_bar(meta,"external_source_RNA_members","external_only_fixed_beta","Fig. S7 | Full predefined bulk-candidate audit","External-only fixed-effect β","FigS7_bulk_full_candidate_audit",COLORS["blue"])
    program=table("104_oep001105_multioics_state_validation","program_score_cross_omic_correlations.tsv"); simple_bar(program.dropna(subset=["rho"]).assign(pair=lambda frame:frame.left_layer+" → "+frame.right_layer),"pair","rho","Fig. S8 | Programme-level cross-omic readout availability","Spearman ρ","FigS8_program_cross_omic_availability",COLORS["protein"])
    clinical=table("104_oep001105_multioics_state_validation","oep001105_vm_clinical_association_tests.tsv"); simple_bar(clinical.assign(test=clinical.test.str.replace("Spearman: ","",regex=False)),"test","fdr_bh_all_prespecified_tests","Fig. S9 | All prespecified surgery-period clinical tests","BH-FDR","FigS9_clinical_covariate_audit",COLORS["grey"])
    genomic=table("105_oep001105_genomic_state_validation","prespecified_genomic_features_protein_directional_check.tsv"); simple_bar(genomic,"feature","RNA_effect","Fig. S10 | Prespecified genomic-feature RNA effects","RNA effect","FigS10_prespecified_genomic_audit",COLORS["orange"])
    response=table("95_gse255058_processed_fpkm_response_analysis","gse255058_processed_fpkm_prespecified_response_test*.tsv"); simple_bar(response,"endpoint","mann_whitney_exact_p","Fig. S11 | GSE255058 treatment-response negative audit","Exact Mann–Whitney P value","FigS11_treatment_response_negative_audit",COLORS["grey"])
    cox=table("100_gse244807_vm_anchor_os_cox_sensitivity","vm_anchor_os_cox_primary_and_clinical_sensitivit*.tsv"); simple_bar(cox,"model","vm_hr_per_1sd","Fig. S12 | GSE244807 OS sensitivity analysis","Hazard ratio per 1-SD score","FigS12_OS_sensitivity_audit",COLORS["grey"])
    spatial=table("68_gse316402_spatial_vm_myeloid_association","myeloid_state_association_summary.tsv"); simple_bar(spatial,"myeloid_state","median_patient_rho","Fig. S13 | BTC spatial extrapolation audit","Median patient-level neighborhood ρ","FigS13_spatial_negative_audit",COLORS["grey"])
    virtual=table("101_cross_cohort_vm_grn_virtual_perturbation","dual_algorithm_virtual_perturbation_target_audit.tsv"); simple_bar(virtual,"TF","minimum_cohort_median_VM_anchor_FDR","Fig. S14 | Dual-algorithm virtual perturbation audit","Minimum cohort median anchor FDR","FigS14_virtual_perturbation_audit",COLORS["grey"])

def main():
    SOURCE.mkdir(parents=True,exist_ok=True); style(); fig5(); fig6(); supplements()

if __name__=="__main__": main()
