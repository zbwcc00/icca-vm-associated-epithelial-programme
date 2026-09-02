# iCCA VM-associated epithelial programme: reproducibility release

This directory is a repository-ready release for the final manuscript assets. It contains the scripts required to regenerate the delivered figures, supplementary Fig. S15 and the assembled Markdown manuscript from the frozen result tables. It does not redistribute public raw or processed data.

## Scope and evidence boundary

The release supports a computational state-map analysis only. It does not establish morphological or perfused vasculogenic mimicry, secretion, ligand–receptor binding, cell–cell communication, causality, treatment resistance, prognosis or a therapeutic dependency.

## Project layout

The scripts expect the following directories at the project root:

```text
analysis/
figures/manuscript_v1/
manuscript_draft/
results/
data/                 # public inputs; excluded from this release
```

Set `ICCA_VM_PROJECT_ROOT` to the directory holding this layout when running the scripts outside the original project location. Without the variable, scripts discover the project root relative to their own `analysis/` directory.

## Environment

Create the pinned Python environment from `environment.lock.yml`, then install the Word-export dependency with `npm install docx@9.6.1`. The exact local environment used Python 3.11.9 and Node.js 24.19.0.

## Rebuild order

1. Acquire the public inputs listed in `SOURCE_DATA_POLICY.md` and recreate the frozen `results/` tables with the full upstream analysis pipeline.
2. From the project root, run `python analysis/130_build_final_submission_assets.py`.
3. Run `node analysis/126_build_submission_docx.js` to regenerate the Word manuscript. The exporter preserves each PNG's native aspect ratio.
4. Check generated files and checksums against `RELEASE_MANIFEST.tsv` before archiving a tagged release.

`analysis/` contains only the asset-generation scripts required for the manuscript package. The historical exploratory scripts in the parent project are retained as provenance and are deliberately not presented as a clean rerunnable pipeline.
