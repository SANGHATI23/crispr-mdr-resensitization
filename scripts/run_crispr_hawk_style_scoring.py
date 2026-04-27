import pandas as pd
from pathlib import Path

inp = Path("results_panstrain/all_panstrain_guide_candidates.csv")
outdir = Path("results_external/crispr_hawk")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inp)

# CRISPR-HAWK-style proxy:
# variant/haplotype-aware stability = pan-strain conservation
# safety = specificity / low off-target burden
# activity = on-target score
# design quality = GC closeness to ideal 50%

df["gc_penalty"] = (df["gc_content"] - 50).abs()
df["gc_quality"] = (100 - df["gc_penalty"]).clip(lower=0)

df["crispr_hawk_proxy_score"] = (
    0.40 * df["conservation_score"] +
    0.30 * df["specificity_score"] +
    0.20 * df["on_target_score"] +
    0.10 * df["gc_quality"]
)

df = df.sort_values("crispr_hawk_proxy_score", ascending=False)

df["rank_crispr_hawk_proxy"] = range(1, len(df) + 1)

cols = [
    "gene", "position", "strand", "spacer", "pam",
    "gc_content", "on_target_score", "specificity_score",
    "conservation_score", "final_score",
    "crispr_hawk_proxy_score", "rank_crispr_hawk_proxy"
]

df[cols].to_csv(outdir / "crispr_hawk_style_output.csv", index=False)

for n in [10, 20, 50]:
    df.head(n)[cols].to_csv(outdir / f"top{n}_crispr_hawk_style.csv", index=False)

# Compare against your original pan-strain final_score ranking
base = pd.read_csv(inp).sort_values("final_score", ascending=False)
hawk = df

rows = []
for n in [10, 20, 50]:
    base_set = set(base.head(n)["spacer"])
    hawk_set = set(hawk.head(n)["spacer"])
    overlap = len(base_set & hawk_set)
    rows.append({
        "top_n": n,
        "overlap_count": overlap,
        "overlap_percent": round(overlap / n * 100, 2)
    })

summary = pd.DataFrame(rows)
summary.to_csv(outdir / "crispr_hawk_overlap_summary.csv", index=False)

print("\nCRISPR-HAWK-style benchmark completed.")
print(summary.to_string(index=False))
print("\nOutputs saved in:", outdir)
