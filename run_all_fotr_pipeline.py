#!/usr/bin/env python3
"""
Run the core FOTR-CRISPR BIBM pipeline.

This script executes the completed analysis modules in order:
1. RNA structural accessibility
2. FOTR v1 scoring
3. CRISOT batch merge
4. CRISOT guide-level summary
5. CRISOT-to-FOTR spacer-level merge
6. Functional context annotation
7. FOTR v2 functional-context scoring
8. FOTR v2 ablation analysis
9. mecA FOTR v2 case study

Important:
- RNAfold/ViennaRNA must be installed and available in the active environment.
- This script assumes it is run from the repository root.
"""

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


STEPS = [
    ("RNA structural accessibility", "analysis/rna_structure/score_rna_structure.py"),
    ("FOTR v1 guide prioritization", "analysis/fotr/score_fotr_v1.py"),
    ("CRISOT batch merge", "analysis/external_merge/merge_crisot_batches.py"),
    ("CRISOT guide-level summary", "analysis/external_merge/summarize_crisot_by_guide.py"),
    ("CRISOT-to-FOTR spacer-level merge", "analysis/external_merge/merge_crisot_into_fotr.py"),
    ("Functional context annotation v1", "analysis/functional_annotation/annotate_functional_context_v1.py"),
    ("FOTR v2 functional-context scoring", "analysis/fotr_v2/score_fotr_v2_functional_context.py"),
    ("FOTR v2 ablation analysis", "analysis/ablation/ablation_fotr_v2.py"),
    ("mecA FOTR v2 case study", "analysis/case_study/case_study_meca_fotr_v2.py"),
]


def run_step(step_name: str, script_path: str) -> None:
    script = ROOT / script_path

    if not script.exists():
        raise FileNotFoundError(f"Missing script: {script_path}")

    print("\n" + "=" * 80)
    print(f"Running: {step_name}")
    print(f"Script:  {script_path}")
    print("=" * 80)

    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=ROOT,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Step failed: {step_name}")


def main() -> None:
    print("Starting FOTR-CRISPR reproducibility pipeline...")
    print(f"Repository root: {ROOT}")

    for step_name, script_path in STEPS:
        run_step(step_name, script_path)

    print("\n" + "=" * 80)
    print("FOTR-CRISPR pipeline completed successfully.")
    print("=" * 80)


if __name__ == "__main__":
    main()