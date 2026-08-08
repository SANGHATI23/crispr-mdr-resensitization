from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = PROJECT_ROOT / "results_external" / "crisot" / "master_batch"
OUTPUT_FILE = (
    PROJECT_ROOT
    / "results_external"
    / "crisot"
    / "merged"
    / "crisot_master_merged.csv"
)


def main() -> None:
    files = sorted(INPUT_DIR.glob("crisot_*.csv"))

    if not files:
        raise FileNotFoundError(f"No CRISOT batch files found in {INPUT_DIR}")

    frames = []

    for file in files:
        df = pd.read_csv(file)
        df["source_file"] = file.name
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)

    merged["sgRNA"] = merged["sgRNA"].astype(str).str.upper().str.strip()
    merged["Target"] = merged["Target"].astype(str).str.upper().str.strip()

    merged["spacer"] = merged["Target"].str[:20]
    merged["pam"] = merged["Target"].str[20:23]
    merged["guide_23mer"] = merged["spacer"] + merged["pam"]

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUTPUT_FILE, index=False)

    print(f"Input files merged: {len(files)}")
    print(f"Total CRISOT rows: {len(merged)}")
    print(f"Unique guide targets: {merged['guide_23mer'].nunique()}")
    print(f"Output saved: {OUTPUT_FILE}")

    print()
    print("First 5 rows:")
    print(
        merged[
            [
                "source_file",
                "spacer",
                "pam",
                "guide_23mer",
                "Mutation",
                "CRISOT-Score",
                "CRISOT-Spec",
                "delta_Spec",
            ]
        ]
        .head()
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()