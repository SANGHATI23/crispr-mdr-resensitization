from pathlib import Path
import re
import subprocess

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results_rna" / "all_guides_full_sgrna.csv"
OUTPUT_FILE = PROJECT_ROOT / "results_rna" / "all_guides_rna_structure.csv"


def run_rnafold(sequence: str) -> tuple[str, float]:
    """
    Run RNAfold on one RNA/DNA sequence.
    RNAfold automatically converts T to U in output.
    Returns dot-bracket structure and MFE.
    """
    result = subprocess.run(
        ["RNAfold", "--noPS"],
        input=sequence,
        text=True,
        capture_output=True,
        check=True,
    )

    lines = result.stdout.strip().splitlines()

    if len(lines) < 2:
        raise ValueError(f"Unexpected RNAfold output:\n{result.stdout}")

    structure_line = lines[1].strip()

    mfe_match = re.search(r"\(([-+]?\d+\.\d+)\)", structure_line)

    if not mfe_match:
        raise ValueError(f"Could not parse MFE from:\n{structure_line}")

    mfe = float(mfe_match.group(1))
    structure = structure_line.split()[0]

    return structure, mfe


def paired_fraction(structure: str, start: int, end: int) -> float:
    """
    Calculate fraction of paired bases in a 0-based interval.
    Paired bases are represented by '(' or ')'.
    """
    region = structure[start:end]

    if len(region) == 0:
        return 0.0

    paired_count = sum(1 for char in region if char in "()")

    return paired_count / len(region)


def longest_hairpin_stem(structure: str, start: int, end: int) -> int:
    """
    Estimate longest continuous paired stretch inside spacer.
    This is a simple warning feature, not experimental proof.
    """
    region = structure[start:end]

    longest = 0
    current = 0

    for char in region:
        if char in "()":
            current += 1
            longest = max(longest, current)
        else:
            current = 0

    return longest


def classify_accessibility(seed_paired_fraction: float, spacer_paired_fraction: float) -> str:
    """
    Conservative structural-accessibility classification.
    Lower paired fraction means more accessible.
    """
    if seed_paired_fraction >= 0.70 or spacer_paired_fraction >= 0.65:
        return "Structurally_Risky"
    elif seed_paired_fraction >= 0.40 or spacer_paired_fraction >= 0.45:
        return "Moderate"
    else:
        return "Accessible"


def main() -> None:
    df = pd.read_csv(INPUT_FILE)

    required_columns = {"gene", "spacer", "full_sgrna_sequence"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    structures = []
    mfe_values = []
    seed_paired_fractions = []
    spacer_paired_fractions = []
    scaffold_paired_fractions = []
    spacer_accessibility_scores = []
    longest_spacer_hairpins = []
    rna_accessibility_classes = []

    for _, row in df.iterrows():
        full_seq = str(row["full_sgrna_sequence"]).upper().strip()

        structure, mfe = run_rnafold(full_seq)

        # sgRNA layout:
        # spacer = first 20 nt
        # seed region = first 10 nt of spacer
        # scaffold = everything after position 20
        seed_pf = paired_fraction(structure, 0, 10)
        spacer_pf = paired_fraction(structure, 0, 20)
        scaffold_pf = paired_fraction(structure, 20, len(structure))

        spacer_accessibility = 1.0 - spacer_pf
        longest_hairpin = longest_hairpin_stem(structure, 0, 20)

        accessibility_class = classify_accessibility(
            seed_paired_fraction=seed_pf,
            spacer_paired_fraction=spacer_pf,
        )

        structures.append(structure)
        mfe_values.append(mfe)
        seed_paired_fractions.append(seed_pf)
        spacer_paired_fractions.append(spacer_pf)
        scaffold_paired_fractions.append(scaffold_pf)
        spacer_accessibility_scores.append(spacer_accessibility)
        longest_spacer_hairpins.append(longest_hairpin)
        rna_accessibility_classes.append(accessibility_class)

    df["rnafold_structure"] = structures
    df["full_sgrna_mfe"] = mfe_values
    df["seed_paired_fraction"] = seed_paired_fractions
    df["spacer_paired_fraction"] = spacer_paired_fractions
    df["scaffold_paired_fraction"] = scaffold_paired_fractions
    df["spacer_accessibility_score"] = spacer_accessibility_scores
    df["longest_spacer_hairpin"] = longest_spacer_hairpins
    df["rna_accessibility_class"] = rna_accessibility_classes

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Input file: {INPUT_FILE}")
    print(f"Output file: {OUTPUT_FILE}")
    print(f"Total guides scored: {len(df)}")
    print("\nRNA accessibility class counts:")
    print(df["rna_accessibility_class"].value_counts())


if __name__ == "__main__":
    main()