import pandas as pd
from pathlib import Path
import subprocess

INPUT = Path("../results_panstrain/all_panstrain_guide_candidates.csv")
OUTPUT = Path("../results_external/cfd/cfd_panel_scored.csv")

BASES = ["A", "C", "G", "T"]


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
    return None


def mutate_base(seq, pos, new_base):
    return seq[:pos] + new_base + seq[pos + 1:]


def build_single_mismatch_panel(spacer, pam="GGG"):
    wt_seq = spacer + pam
    panel = []

    # mutate only guide region: first 20 nt
    for pos in range(20):
        original = spacer[pos]
        for base in BASES:
            if base != original:
                mutated_spacer = mutate_base(spacer, pos, base)
                off_seq = mutated_spacer + pam
                panel.append({
                    "position": pos + 1,
                    "ref_base": original,
                    "alt_base": base,
                    "off_seq": off_seq
                })
    return wt_seq, panel


def summarize_scores(scores):
    valid = [s for s in scores if s is not None]
    if not valid:
        return None, None, None
    return min(valid), sum(valid) / len(valid), max(valid)


def main():
    df = pd.read_csv(INPUT)

    min_scores = []
    mean_scores = []
    max_scores = []
    panel_counts = []

    for i, row in df.iterrows():
        spacer = row["spacer"]
        wt_seq, panel = build_single_mismatch_panel(spacer)

        scores = []
        for item in panel:
            score = run_cfd(wt_seq, item["off_seq"])
            scores.append(score)

        min_cfd, mean_cfd, max_cfd = summarize_scores(scores)

        min_scores.append(min_cfd)
        mean_scores.append(mean_cfd)
        max_scores.append(max_cfd)
        panel_counts.append(len(panel))

        if i % 10 == 0:
            print(f"Processed {i} guides")

    df["cfd_panel_min"] = min_scores
    df["cfd_panel_mean"] = mean_scores
    df["cfd_panel_max"] = max_scores
    df["cfd_panel_n"] = panel_counts

    df["rank_cfd_panel_mean"] = df["cfd_panel_mean"].rank(
        ascending=False, method="min"
    )
    df["rank_cfd_panel_min"] = df["cfd_panel_min"].rank(
        ascending=False, method="min"
    )
    df["rank_final"] = df["final_score"].rank(
        ascending=False, method="min"
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT, index=False)

    print("Done. Saved to:", OUTPUT)
    print(df[[
        "gene", "spacer", "final_score",
        "cfd_panel_min", "cfd_panel_mean", "cfd_panel_max",
        "rank_final", "rank_cfd_panel_mean"
    ]].sort_values("rank_final").head(10).to_string(index=False))


if __name__ == "__main__":
    main()