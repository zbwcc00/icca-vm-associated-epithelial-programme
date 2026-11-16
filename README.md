# iCCA VM-related epithelial programme: corrected revision release

This staging release corresponds to the corrected resubmission package for *Cross-cohort definition and multi-omic characterization of a six-gene VM-related epithelial transcriptional programme in intrahepatic cholangiocarcinoma*.

## Corrected truth policy

The complete six-gene raw-UMI reconstruction is the sole analytical truth for this corrected revision. All label-dependent scores, q80/q90 labels, pseudobulk contrasts, SecAct funnels, bulk candidate audits, figures and supplementary tables were regenerated or withdrawn from the public release namespace accordingly.

Corrected SecAct funnels are 657 → 268 → 58 → 26 at q80 and 657 → 244 → 56 → 26 at q90. Corrected bulk audit: 24 groups evaluable in all three cohorts, 12 positive in all three and 17 with the same direction. These are association and calibration analyses, not evidence of morphology, perfusion, secretion, receptor engagement, causal communication, prognosis, treatment resistance or a therapeutic target.

## Contents

- `manuscript/` contains the corrected clean manuscript. Peer-review correspondence and the red-marked manuscript are retained only in the journal submission package.
- `figures/` contains revised main and supplementary figure assets.
- `source_data/` contains corrected panel-level source tables; superseded tables are not exposed here.
- `supplementary_tables/S7_corrected/` contains corrected Reviewer 3 audit tables.
- `analysis/reviewer3_corrected_truth/` contains correction-audit provenance scripts.

## Reproduction boundary

This is an auditable manuscript-asset and corrected-source-data release. Public raw and processed matrices are not redistributed. Accessions are listed in `metadata/data_accessions.tsv`. Some provenance scripts require the frozen upstream results directory and are not a one-command raw-data reanalysis.

## Archived release

This corrected analysis release is v1.1.1. Cite the all-version Zenodo concept DOI [10.5281/zenodo.22255463](https://doi.org/10.5281/zenodo.22255463); it resolves to the latest immutable release. The corresponding GitHub release is [v1.1.1](https://github.com/zbwcc00/icca-vm-associated-epithelial-programme/releases/tag/v1.1.1).
