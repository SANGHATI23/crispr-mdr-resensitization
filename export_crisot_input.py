import os
import pandas as pd

INPUT_CSV = "results/top20_global_guides.csv"
OUTPUT_TXT = "results_external/crisot/crisot_input_sgrnas.txt"

def main():
    df = pd.read_csv(INPUT_CSV)

    required = {"gene", "spacer", "pam"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df["Target sequence"] = (
        df["spacer"].astype(str).str.upper()
        + df["pam"].astype(str).str.upper()
    )

    os.makedirs("results_external/crisot", exist_ok=True)

    out = df[["gene", "spacer", "pam", "Target sequence"]].copy()
    out.to_csv(OUTPUT_TXT, sep="\t", index=False)

    print("Saved:", OUTPUT_TXT)
    print("Rows:", len(out))
    print(out.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
