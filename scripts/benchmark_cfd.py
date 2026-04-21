import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

INPUT_CSV = REPO_ROOT / "results_panstrain" / "all_panstrain_guide_candidates.csv"
OUTPUT_CSV = REPO_ROOT / "results_external" / "cfd" / "cfd_scores.csv"


def simple_cfd_score(spacer: str, off_target: str) -> float:
    """
    Simplified CFD-like mismatch penalty
    """
    score = 1.0
    for a, b in zip(spacer, off_target):
        if a != b:
            score *= 0.8  # penalty per mismatch
    return score


def main():
    df = pd.read_csv(INPUT_CSV)

    if "spacer" not in df.columns:
        raise ValueError("Missing spacer column")

    # For now: simulate off-targets as mutated versions
    df["mock_offtarget"] = df["spacer"].apply(lambda x: x[:-1] + "A")

    df["cfd_score"] = df.apply(
        lambda row: simple_cfd_score(row["spacer"], row["mock_offtarget"]),
        axis=1
    )

    df["cfd_score_100"] = (df["cfd_score"] * 100).round(2)

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)

    print("CFD scoring complete")
    print(df[["spacer", "cfd_score_100"]].head())


if __name__ == "__main__":
    main()
