import os
from pathlib import Path

ROOT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
DRAFT = ROOT / "manuscript_draft"
OUT = DRAFT / "09_full_manuscript_draft_v1.md"


def section_body(filename: str, drop_editorial_notes: bool = False) -> str:
    text = (DRAFT / filename).read_text(encoding="utf-8")
    lines = text.splitlines()
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
    text = "\n".join(lines).strip()
    if drop_editorial_notes and "## Editorial supervision notes" in text:
        text = text.split("## Editorial supervision notes", 1)[0].rstrip()
    return text


def final_references() -> str:
    """Remove the unrelated legacy entry and normalise final Vancouver punctuation."""
    references = (ROOT / "results" / "111_reference_search" / "final_verified_references_vancouver.txt").read_text(encoding="utf-8").splitlines()
    kept = [line for line in references if not line.startswith("24. Cong Y,")]
    output = []
    for number, line in enumerate(kept, start=1):
        _, body = line.split(". ", 1)
        output.append(f"{number}. {body.replace('..', '.')}")
    return "\n".join(output)


title = "A patient-aware multi-omic analysis of a vasculogenic-mimicry-associated epithelial transcriptional programme in intrahepatic cholangiocarcinoma"

parts = [
    f"# {title}",
    "",
    "**Authors:** Bowen Zheng, Guanghua Xie, Depei Kong, Jianlei Chen, Ye Wang, Hao Li*",
    "",
    "**Affiliation:** Division of Hepatobiliary Pancreatic Surgery, The Affiliated Hospital of Yanbian University, Yanji 133000, China",
    "",
    "**Author emails and ORCID iDs:** Bowen Zheng (13803722355@163.com; ORCID: 0009-0006-8304-7905); Guanghua Xie (ghxie@ybu.edu.cn; ORCID: 0000-0002-6461-6707); Depei Kong (18764537280@163.com; ORCID: 0009-0008-1904-7198); Jianlei Chen (chenjl421@163.com; ORCID: 0009-0007-5316-856X); Ye Wang (mercy9909@163.com; ORCID: 0009-0000-1282-6123); Hao Li (lih@ybu.edu.cn; ORCID: 0009-0007-7394-1273).",
    "",
    "*Corresponding author: Hao Li (lih@ybu.edu.cn).",
    "",
    "**Manuscript type:** Research article (computational multi-omic analysis)",
    "",
    "**Keywords:** intrahepatic cholangiocarcinoma; vasculogenic-mimicry-associated state; single-cell RNA sequencing; proteogenomics; tumour microenvironment; somatic alteration",
    "",
    "## Abstract",
    section_body("08_abstract_draft.md"),
    "",
    "## Introduction",
    section_body("07_introduction_draft.md"),
    "",
    "## Methods",
    section_body("02_methods_draft.md"),
    "",
    "## Results",
    section_body("03_results_draft.md"),
    "",
    "## Discussion",
    section_body("06_discussion_draft.md", drop_editorial_notes=True),
    "",
    "## Data availability",
    "All analyses use de-identified public resources. Accession-level source datasets, cohort eligibility records and sample-selection tables are listed in Supplementary Tables S1–S3 and will be deposited with the analysis repository. Public raw and processed data remain available from GEO, BioStudies, BioSino/NODE accession OEP001105 and the originating repositories under their respective access conditions. This manuscript is not submission-ready until a public repository URL and archived release DOI are inserted here.",
    "",
    "## Code availability",
    "The repository-ready manuscript-asset release contains source tables, execution manifests, fixed random seeds and a pinned environment specification, with relative project-root discovery or a documented environment-variable override. A versioned public repository and tagged archival release will be created before submission. This manuscript is not submission-ready until the repository URL, commit/tag and release DOI are inserted here.",
    "",
    "## Ethics statement",
    "This study used de-identified public data generated under the ethics approvals of the originating studies. No new human samples, interventions or identifiable clinical data were collected by the authors.",
    "",
    "## Author contributions",
    "Bowen Zheng: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Visualization, Writing – original draft. Guanghua Xie: Methodology, Data curation, Validation, Resources, Writing – review & editing. Depei Kong: Data curation, Validation, Investigation, Resources, Writing – review & editing. Jianlei Chen: Software, Formal analysis, Visualization, Data curation. Ye Wang: Resources, Data curation, Validation, Writing – review & editing. Hao Li: Supervision, Project administration, Writing – review & editing, Final approval of the manuscript. All authors read and approved the final manuscript.",
    "",
    "## Funding",
    "This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.",
    "",
    "## Competing interests",
    "The authors declare no competing interests.",
    "",
    "## Acknowledgements",
    "None.",
    "",
    "## AI assistance disclosure",
    "Language and structure were assisted by an AI system under human author supervision. All analyses, source data, statistical results, interpretations and final claims require author verification before submission. The AI system is not an author.",
    "",
    "## Figure legends",
    section_body("01_figure_legends_draft.md"),
    "",
    "## Supplementary materials guide",
    "Supplementary Figures S1–S15 provide quality-control, scoring, source-RNA, bulk, multi-omic, genomic, treatment-response, survival, spatial, robustness and virtual-perturbation audits. Their exact panel descriptions and limitations are included in the Figure legends section above. Supplementary source tables and the claim-level evidence audit will be included in the versioned public repository release; repository URL and release DOI remain to be inserted before submission.",
    "",
    "## References",
    final_references(),
]

OUT.write_text("\n".join(parts).strip() + "\n", encoding="utf-8")
print(OUT)
