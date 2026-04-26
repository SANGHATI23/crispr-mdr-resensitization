import os
import pandas as pd

INPUT_CSV = "results/top20_global_guides.csv"
OUTPUT_CSV = "results_external/azimuth/azimuth_input.csv"

def build_30mer(spacer, pam):
    # 4 nt upstream + spacer + PAM + 3 nt downstream
    return "AAAA" + spacer + pam + "GGG"

def main():
    df = pd.read_csv(INPUT_CSV)

    df["context_30mer"] = df.apply(
        lambda r: build_30mer(r["spacer"], r["pam"]),
        axis=1
    )

    os.makedirs("results_external/azimuth", exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)

    print("Saved:", OUTPUT_CSV)
    print("Rows:", len(df))
    print(df[["gene","spacer","context_30mer"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
