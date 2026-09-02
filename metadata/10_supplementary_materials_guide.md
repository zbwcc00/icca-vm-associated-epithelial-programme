# Supplementary materials guide

## Supplementary figures

| Item | Content | Evidence role | Boundary |
|---|---|---|---|
| Fig. S1 | Single-cell library QC and `scDblFinder` doublet audit | provenance and QC | descriptive; not patient inference |
| Fig. S2 | CNV-reference sensitivity and malignant-epithelial annotation | classification robustness | operational label, not pathology gold standard |
| Fig. S3 | UCell/AUCell/module-score concordance and consensus frequencies | score reliability | descriptive frequencies |
| Fig. S4 | Myeloid and stromal state coverage | eligibility for source analysis | low-coverage/mixed states excluded |
| Fig. S5 | External patient-by-state SecAct source-RNA audit | qualification of 66 candidates | two-patient compatibility audit; no P value |
| Fig. S6 | GSE201425 external myeloid identity transfer | external identity support | three states supported; FCGR3A exclusion-only |
| Fig. S7 | Full bulk candidate and meta-analysis audit | transparency of candidate triage | group 014 direction reversal retained |
| Fig. S8 | Cross-omic availability and programme correlations | missingness and concordance | no complete phospho-programme |
| Fig. S9 | Clinical covariate and HC3 sensitivity audit | confounding boundary | not clinical independence or prognosis |
| Fig. S10 | Prespecified mutation/CNV audit | genomic transparency | RNA event set n=249; matched protein directional checks n=204; genome-wide CNV scan negative |
| Fig. S11 | GSE255058 treatment-response audit | negative endpoint | 9 responders versus 9 non-responders; underpowered |
| Fig. S12 | GSE244807 OS and specimen-type sensitivity | outcome boundary | unadjusted association not robust to tissue adjustment |
| Fig. S13 | GSE316402 BTC spatial extrapolation audit | transferability boundary | BTC, not iCCA; corrected sign tests negative |
| Fig. S14 | GRN, dual virtual perturbation, PRISM and CRISPR audit | non-actionability boundary | no validated upstream target or clinical drug claim |
| Fig. S15 | SecAct global-FDR, published-programme, fixed-seed expression-matched random-set, anchor-ablation and leave-one-cohort-out robustness | calibration and sensitivity of molecular coordinate | does not establish VM specificity, morphology, causality or independent clinical validation |

## Supplementary source data

The versioned public repository contains the numerical source tables used for the displayed panels, final figures, supplementary figures, supplementary tables, analysis provenance scripts, an execution manifest and the environment lockfile. It does not redistribute raw or processed source matrices. Source tables preserve patient/sample identifiers in the form supplied by the public repositories and include the analysis unit, sample count and correction family wherever available. Public raw and processed source data remain available from their originating repositories under the access conditions stated by those repositories.

## Supplementary reporting notes

1. Single-cell coverage summaries use patients; receiver-side contrasts use paired patient-or-library pseudobulks. Libraries are not interpreted as independent patients, and cell numbers shown in plots are descriptive.
2. SecAct group counts (657, 249, 66 and 30) are nested filtering stages, not independent biological replicates.
3. OEP001105 RNA, protein, phosphopeptide and genomic layers are matched within one project and are not independent cohort validation.
4. The GSE244807 OS analysis is retained as a sensitivity audit because the estimate changes materially with specimen type.
5. Negative treatment, spatial, virtual-perturbation, CRISPR and pharmacogenomic analyses are not omitted from the evidence record.
6. Supplementary Table S1 lists single-cell library quality-control thresholds and doublet summaries; Supplementary Table S2 lists cohort, patient/library/cell counts, operational-state criteria and pseudobulk eligibility; Supplementary Table S3 lists SecAct testing families, source-RNA rules, multi-omic matching and missing-data handling; Supplementary Table S4 contains the Fig. S15 calibration source data, including patient-level anchor-ablation effects.
