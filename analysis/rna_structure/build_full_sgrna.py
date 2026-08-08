from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "results_panstrain"
    / "all_panstrain_guide_candidates.csv"
)

OUTPUT_DIR = PROJECT_ROOT / "results_rna"
OUTPUT_FILE = OUTPUT_DIR / "all_guides_full_sgrna.csv"

SPCAS9_SGRNA_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG"
    "AAAAAGTGGCACCGAGTCGGTGCTTTTT"
)


def validate_spacers(df: pd.DataFrame) -> None:
    required_columns = {"gene", "spacer", "pam", "final_score"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}"
        )

    df["spacer"] = df["spacer"].astype(str).str.upper().str.strip()

    invalid_lengths = df[df["spacer"].str.len() != 20]

    if not invalid_lengths.empty:
        raise ValueError(
            f"Found {len(invalid_lengths)} spacers that are not 20 nt long."
        )

    invalid_sequences = df[
        ~df["spacer"].str.fullmatch(r"[ACGT]{20}")
    ]

    if not invalid_sequences.empty:
        raise ValueError(
            f"Found {len(invalid_sequences)} spacers with invalid characters."
        )


def main() -> None:
    df = pd.read_csv(INPUT_FILE)

    validate_spacers(df)

    df["full_sgrna_sequence"] = (
        df["spacer"] + SPCAS9_SGRNA_SCAFFOLD
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Input file: {INPUT_FILE}")
    print(f"Total guides: {len(df)}")
    print(f"Spacer length: {df['spacer'].str.len().iloc[0]} nt")
    print(
        "Scaffold length: "
        f"{len(SPCAS9_SGRNA_SCAFFOLD)} nt"
    )
    print(
        "Full sgRNA length: "
        f"{df['full_sgrna_sequence'].str.len().iloc[0]} nt"
    )
    print(f"Output saved to: {OUTPUT_FILE}")

    print("\nFirst full sgRNA:")
    print(df.loc[0, "full_sgrna_sequence"])


if __name__ == "__main__":
    main()