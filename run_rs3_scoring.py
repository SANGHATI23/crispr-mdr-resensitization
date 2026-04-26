import os
import pandas as pd
from rs3.seq import predict_seq

INPUT_CSV = "results_external/azimuth/azimuth_input.csv"
OUTPUT_CSV = "results_external/rs3/rs3_scored.csv"

def main():
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    if "context_30mer" not in df.columns:
        raise ValueError("Column 'context_30mer' not found in input CSV")

    seqs = df["context_30mer"].astype(str).tolist()

    print("Loaded rows:", len(df))
    print("Scoring with Rule Set 3...")

    preds = predict_seq(
        seqs,
        sequence_tracr="Hsu2013"
    )

    df["rs3_score"] = preds
    df = df.sort_values("rs3_score", ascending=False).reset_index(drop=True)
    df["rank_rs3"] = df.index + 1

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)

    print("\n==============================")
    print("RS3 SCORING COMPLETE")
    print("==============================")
    print("Output file:", OUTPUT_CSV)
    print("\nTop 10 rows:")
    print(df.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
