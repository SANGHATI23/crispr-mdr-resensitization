import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

OUTDIR = Path("results/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv("results/final_comparison_table.csv")

# Auto-detect columns
guide_col = "guide_sequence" if "guide_sequence" in df.columns else "guide"
final_col = "final_score"
mit_col = "MIT_score" if "MIT_score" in df.columns else "mit_score"
crisot_col = "crisot_score"
cfd_col = "CFD_score" if "CFD_score" in df.columns else "cfd_score"

# Keep only rows with needed scores
plot_df = df[[guide_col, final_col, mit_col, crisot_col, cfd_col]].dropna()

# Normalize MIT if needed
# MIT may already be 0-1, while others may be 0-100
if plot_df[mit_col].max() <= 1.5:
    plot_df["MIT_score_scaled"] = plot_df[mit_col] * 100
else:
    plot_df["MIT_score_scaled"] = plot_df[mit_col]

# Figure 1: Composite vs CRISOT
plt.figure(figsize=(7, 5))
plt.scatter(plot_df[final_col], plot_df[crisot_col])
plt.xlabel("Composite Score")
plt.ylabel("CRISOT Score")
plt.title("Composite Score vs CRISOT Score")
plt.tight_layout()
plt.savefig(OUTDIR / "composite_vs_crisot.png", dpi=300)
plt.close()

# Figure 2: Composite vs MIT
plt.figure(figsize=(7, 5))
plt.scatter(plot_df[final_col], plot_df["MIT_score_scaled"])
plt.xlabel("Composite Score")
plt.ylabel("MIT Score (scaled 0-100)")
plt.title("Composite Score vs MIT Score")
plt.tight_layout()
plt.savefig(OUTDIR / "composite_vs_mit.png", dpi=300)
plt.close()

# Figure 3: Composite vs CFD
plt.figure(figsize=(7, 5))
plt.scatter(plot_df[final_col], plot_df[cfd_col])
plt.xlabel("Composite Score")
plt.ylabel("CFD Score")
plt.title("Composite Score vs CFD Score")
plt.tight_layout()
plt.savefig(OUTDIR / "composite_vs_cfd.png", dpi=300)
plt.close()

# Figure 4: Top 15 composite guide comparison
top = plot_df.sort_values(final_col, ascending=False).head(15).copy()
top["short_guide"] = top[guide_col].str[:6] + "..." + top[guide_col].str[-4:]

x = range(len(top))

plt.figure(figsize=(12, 6))
plt.plot(x, top[final_col], marker="o", label="Composite")
plt.plot(x, top["MIT_score_scaled"], marker="o", label="MIT scaled")
plt.plot(x, top[crisot_col] * 100 if top[crisot_col].max() <= 1.5 else top[crisot_col], marker="o", label="CRISOT scaled")
plt.plot(x, top[cfd_col], marker="o", label="CFD")
plt.xticks(x, top["short_guide"], rotation=45, ha="right")
plt.ylabel("Score")
plt.xlabel("Top Composite Guides")
plt.title("Top Composite Guides Across Scoring Systems")
plt.legend()
plt.tight_layout()
plt.savefig(OUTDIR / "top15_guides_cross_model_scores.png", dpi=300)
plt.close()

print("Saved figures to:", OUTDIR)
print("Files:")
for f in OUTDIR.glob("*.png"):
    print("-", f)
