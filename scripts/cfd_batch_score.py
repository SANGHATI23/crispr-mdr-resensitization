import pandas as pd
from pathlib import Path
import subprocess

INPUT = Path("../results_panstrain/all_panstrain_guide_candidates.csv")
OUTPUT = Path("../results_external/cfd/cfd_scored.csv")

def run_cfd(wt, off):
    cmd = [
        "python",
        "cfd-score-calculator.py",
        "--wt", wt,
        "--off", off
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if "CFD score" in result.stdout:
        return float(result.stdout.strip().split(": ")[1])
    else:
        return None

def main():
    df = pd.read_csv(INPUT)

    # Create mock off-target (1 mismatch inside first 20 nt)
    df["off_target"] = df["spacer"].apply(lambda x: x[:18] + "A" + x[19:] + "GGG")

    # Create WT as spacer + PAM (GGG)
    df["wt_seq"] = df["spacer"] + "GGG"

    scores = []
    for i, row in df.iterrows():
        score = run_cfd(row["wt_seq"], row["off_target"])
        scores.append(score)

        if i % 50 == 0:
            print(f"Processed {i}")

    df["cfd_score"] = scores

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT, index=False)

    print("Done. Saved to:", OUTPUT)

if __name__ == "__main__":
    main()