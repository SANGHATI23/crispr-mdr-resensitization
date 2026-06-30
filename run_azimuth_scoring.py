import os
import pandas as pd
import azimuth.model_comparison

INPUT_CSV = "results_external/azimuth/azimuth_input.csv"
OUTPUT_CSV = "results_external/azimuth/azimuth_scored.csv"

def main():
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    if "context_30mer" not in df.columns:
        raise ValueError("Column 'context_30mer' not found in input CSV")

    seqs = df["context_30mer"].astype(str).tolist()

    print("Loaded rows:", len(df))
    print("Scoring with Azimuth...")

    predictions = azimuth.model_comparison.predict(
        seq=seqs,
        aa_cut=None,
        percent_peptide=None
    )

    df["azimuth_score"] = predictions
    df = df.sort_values("azimuth_score", ascending=False).reset_index(drop=True)
    df["rank_azimuth"] = df.index + 1

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)

    print("\n==============================")
    print("AZIMUTH SCORING COMPLETE")
    print("==============================")
    print("Output file:", OUTPUT_CSV)
    print("\nTop 10 rows:")
    print(df.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
