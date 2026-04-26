import pandas as pd
from pathlib import Path

INPUT = Path("results_panstrain/all_panstrain_guide_candidates.csv")
OUTDIR = Path("results_external")
OUTDIR.mkdir(exist_ok=True)

OUTPUT = OUTDIR / "sscrispr_proxy_scores.csv"

def gc_score(gc):
    if 40 <= gc <= 60:
        return 100
    elif 35 <= gc < 40 or 60 < gc <= 70:
        return 80
    else:
        return 50

def homopolymer_penalty(seq):
    penalties = 0
    for base in "ATGC":
        if base * 4 in seq:
            penalties += 10
        if base * 5 in seq:
            penalties += 20
    return penalties

def seed_specificity_score(seq):
    seed = seq[-10:]
    seed_gc = (seed.count("G") + seed.count("C")) / len(seed) * 100

    score = 100
    if seed_gc < 30 or seed_gc > 80:
        score -= 15
    if "TTTT" in seq:
        score -= 20
    score -= homopolymer_penalty(seq)

    return max(0, score)

def sscrispr_proxy_score(row):
    spacer = str(row["spacer"]).upper()
    gc = float(row["gc_content"])

    gc_component = gc_score(gc)
    seed_component = seed_specificity_score(spacer)
    conservation = float(row.get("conservation_score", 100))
    specificity = float(row.get("specificity_score", 100))
    on_target = float(row.get("on_target_score", 0))

    final = (
        0.30 * on_target +
        0.25 * specificity +
        0.20 * conservation +
        0.15 * gc_component +
        0.10 * seed_component
    )

    return round(final, 2)

df = pd.read_csv(INPUT)

df["sscrispr_proxy_score"] = df.apply(sscrispr_proxy_score, axis=1)

df["sscrispr_proxy_rank"] = (
    df["sscrispr_proxy_score"]
    .rank(method="dense", ascending=False)
    .astype(int)
)

df["sscrispr_proxy_class"] = pd.cut(
    df["sscrispr_proxy_score"],
    bins=[-1, 50, 70, 85, 100],
    labels=["Poor", "Moderate", "Good", "Excellent"]
)

df = df.sort_values(["sscrispr_proxy_score", "final_score"], ascending=False)

df.to_csv(OUTPUT, index=False)

print(f"Saved: {OUTPUT}")
print(
    df[
        [
            "gene",
            "spacer",
            "pam",
            "sscrispr_proxy_score",
            "sscrispr_proxy_rank",
            "sscrispr_proxy_class",
        ]
    ].head(20)
)
