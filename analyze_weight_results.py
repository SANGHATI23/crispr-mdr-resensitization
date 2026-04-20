import pandas as pd
from pathlib import Path

# ----------------------------
# Load files
# ----------------------------
w2 = pd.read_csv("weight_sensitivity_results/phase2_weight_summary_all_variations.csv")
g2 = pd.read_csv("weight_sensitivity_results/phase2_per_gene_summary_all_variations.csv")
a2 = pd.read_csv("weight_sensitivity_results/phase2_all_weighted_guides_all_variations.csv")

# Output folder
outdir = Path("results_weight_analysis")
outdir.mkdir(exist_ok=True)

print("\nLoaded files successfully.")
print(f"w2 shape: {w2.shape}")
print(f"g2 shape: {g2.shape}")
print(f"a2 shape: {a2.shape}")

# ----------------------------
# 1. Compare 40/30/30 with nearby weights
# ----------------------------
focus_weights = ["40/30/30", "50/30/20", "30/40/30", "40/40/20", "30/30/40"]
focus = w2[w2["weight_scheme"].isin(focus_weights)].copy()

print("\n=== 40/30/30 AND NEARBY WEIGHTS ===")
print(focus.sort_values("weight_scheme"))

focus.sort_values("weight_scheme").to_csv(
    outdir / "phase2_focus_weights_comparison.csv", index=False
)

# ----------------------------
# 2. Best overall schemes by different criteria
# ----------------------------
best_mean = w2.sort_values("mean_final_score", ascending=False).head(10)
best_excellent = w2.sort_values("excellent_count", ascending=False).head(10)
best_overlap = w2.sort_values("top20_jaccard_vs_40_30_30", ascending=False).head(10)

print("\n=== TOP 10 BY MEAN FINAL SCORE ===")
print(best_mean[[
    "weight_scheme", "mean_final_score", "best_final_score",
    "excellent_count", "top20_jaccard_vs_40_30_30",
    "top_guide_gene", "top_guide_spacer", "top_guide_score"
]])

print("\n=== TOP 10 BY EXCELLENT COUNT ===")
print(best_excellent[[
    "weight_scheme", "excellent_count", "mean_final_score",
    "best_final_score", "top20_jaccard_vs_40_30_30",
    "top_guide_gene", "top_guide_spacer", "top_guide_score"
]])

print("\n=== TOP 10 BY TOP-20 JACCARD VS 40/30/30 ===")
print(best_overlap[[
    "weight_scheme", "top20_overlap_count_vs_40_30_30",
    "top20_jaccard_vs_40_30_30", "mean_final_score",
    "excellent_count", "top_guide_gene", "top_guide_spacer"
]])

best_mean.to_csv(outdir / "phase2_top10_by_mean_score.csv", index=False)
best_excellent.to_csv(outdir / "phase2_top10_by_excellent_count.csv", index=False)
best_overlap.to_csv(outdir / "phase2_top10_by_jaccard.csv", index=False)

# ----------------------------
# 3. Top guide stability across all 66 combinations
# ----------------------------
top_guide_stability = (
    w2.groupby(["top_guide_gene", "top_guide_spacer"])
      .size()
      .reset_index(name="number_of_weight_schemes")
      .sort_values("number_of_weight_schemes", ascending=False)
)

print("\n=== TOP GUIDE STABILITY ACROSS ALL WEIGHTS ===")
print(top_guide_stability)

top_guide_stability.to_csv(outdir / "phase2_top_guide_stability.csv", index=False)

# ----------------------------
# 4. Top 1 guide per gene per weight
# ----------------------------
a2_sorted = a2.sort_values(
    ["weight_scheme", "gene", "final_score_new"],
    ascending=[True, True, False]
).copy()

top1_per_gene = (
    a2_sorted.groupby(["weight_scheme", "gene"])
    .head(1)
    .copy()
)

print("\n=== TOP 1 GUIDE PER GENE PER WEIGHT (first 20 rows) ===")
print(top1_per_gene[[
    "weight_scheme", "gene", "spacer", "final_score_new", "classification_new"
]].head(20))

top1_per_gene.to_csv(outdir / "phase2_top1_per_gene_per_weight.csv", index=False)

# ----------------------------
# 5. For each gene, which spacer becomes top most often?
# ----------------------------
gene_top_stability = (
    top1_per_gene.groupby(["gene", "spacer"])
    .size()
    .reset_index(name="times_ranked_top1")
    .sort_values(["gene", "times_ranked_top1"], ascending=[True, False])
)

print("\n=== PER-GENE TOP GUIDE STABILITY ===")
print(gene_top_stability.head(30))

gene_top_stability.to_csv(outdir / "phase2_per_gene_top1_stability.csv", index=False)

# ----------------------------
# 6. Show the exact row for 40/30/30
# ----------------------------
ref = w2[w2["weight_scheme"] == "40/30/30"].copy()

print("\n=== EXACT 40/30/30 ROW ===")
print(ref)

ref.to_csv(outdir / "phase2_reference_40_30_30.csv", index=False)

print("\nDone. Files written to:")
print(outdir.resolve())