"""Rebuild the final manuscript figures and assembled Markdown from frozen results.

Set ICCA_VM_PROJECT_ROOT to run this driver from another project location.
The Word export remains a separate Node.js step because it uses docx.
"""
from __future__ import annotations

import os
import runpy
import subprocess
import sys
from pathlib import Path


ROOT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
ANALYSIS = ROOT / "analysis"


def run_python(script_name: str) -> None:
    subprocess.run([sys.executable, str(ANALYSIS / script_name)], check=True, cwd=ROOT)


def main() -> None:
    first_figure = runpy.run_path(str(ANALYSIS / "114_build_manuscript_figures_v1.py"))
    first_figure["setup_style"]()
    first_figure["figure_1"]()
    for script in [
        "115_rebuild_manuscript_figures_2_4.py",
        "116_rebuild_figures_5_6_and_supplementary.py",
        "129_build_supplementary_figure_s15.py",
        "117_assemble_full_manuscript.py",
    ]:
        run_python(script)


if __name__ == "__main__":
    main()
