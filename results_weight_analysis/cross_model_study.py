import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ===== Load dataset =====
file_path = REPO_ROOT / "weight_sensitivity_results" / "phase2_all_weighted_guides_all_variations.csv"
df = pd.read_csv(file_path)

# ===== Filter reference weight =====
ref_weight = "40/30/30"
df_ref = df[df["weight_scheme"] == ref_weight].copy()

print("Loaded:", file_path.name)
print("Reference weight:", ref_weight)
print("Filtered shape:", df_ref.shape)

# ===== Create ranks =====
df_ref["rank_on_target"] = df_ref["on_target_score"].rank(ascending=False, method="min")
df_ref["rank_off_target"] = df_ref["specificity_score"].rank(ascending=False, method="min")
df_ref["rank_conservation"] = df_ref["conservation_score"].rank(ascending=False, method="min")
df_ref["rank_final"] = df_ref["final_score_new"].rank(ascending=False, method="min")

# ===== Consensus model =====
df_ref["consensus_score"] = (
    df_ref["on_target_score"] +
    df_ref["specificity_score"] +
    df_ref["conservation_score"]
) / 3.0
df_ref["rank_consensus"] = df_ref["consensus_score"].rank(ascending=False, method="min")

# ===== Save full ranked reference dataset =====
full_ranked_output = RESULTS_DIR / "cross_model_reference_ranked_guides.csv"
df_ref.to_csv(full_ranked_output, index=False)

# ===== Top 10 by each model =====
top_n = 10

top_on_target_df = df_ref.nsmallest(top_n, "rank_on_target").copy()
top_off_target_df = df_ref.nsmallest(top_n, "rank_off_target").copy()
top_conservation_df = df_ref.nsmallest(top_n, "rank_conservation").copy()
top_consensus_df = df_ref.nsmallest(top_n, "rank_consensus").copy()
top_final_df = df_ref.nsmallest(top_n, "rank_final").copy()

# ===== Save top 10 files separately =====
top_on_target_df.to_csv(RESULTS_DIR / "top10_on_target_model.csv", index=False)
top_off_target_df.to_csv(RESULTS_DIR / "top10_off_target_model.csv", index=False)
top_conservation_df.to_csv(RESULTS_DIR / "top10_conservation_model.csv", index=False)
top_consensus_df.to_csv(RESULTS_DIR / "top10_consensus_model.csv", index=False)
top_final_df.to_csv(RESULTS_DIR / "top10_final_composite_model.csv", index=False)

# ===== Combined comparison table using union of all top-10 guides =====
all_top_guides = set(top_on_target_df["guide_key"]) \
    | set(top_off_target_df["guide_key"]) \
    | set(top_conservation_df["guide_key"]) \
    | set(top_consensus_df["guide_key"]) \
    | set(top_final_df["guide_key"])

comparison_cols = [
    "gene",
    "spacer",
    "guide_key",
    "on_target_score",
    "specificity_score",
    "conservation_score",
    "consensus_score",
    "final_score_new",
    "rank_on_target",
    "rank_off_target",
    "rank_conservation",
    "rank_consensus",
    "rank_final",
]

comparison_df = df_ref[df_ref["guide_key"].isin(all_top_guides)][comparison_cols].copy()
comparison_df = comparison_df.sort_values(["rank_final", "rank_consensus", "rank_on_target"])

comparison_output = RESULTS_DIR / "cross_model_top_guides_comparison_table.csv"
comparison_df.to_csv(comparison_output, index=False)

# ===== Disagreement analysis: high on-target but dropped by final =====
high_on_target = df_ref.nsmallest(20, "rank_on_target").copy()
dropped_guides = high_on_target[
    ~high_on_target["guide_key"].isin(top_final_df["guide_key"])
].copy()

dropped_cols = [
    "gene",
    "spacer",
    "guide_key",
    "on_target_score",
    "specificity_score",
    "conservation_score",
    "final_score_new",
    "rank_on_target",
    "rank_final",
]

dropped_guides = dropped_guides[dropped_cols].sort_values(["rank_on_target", "rank_final"])
dropped_output = RESULTS_DIR / "high_on_target_but_dropped_by_final.csv"
dropped_guides.to_csv(dropped_output, index=False)

# ===== Overlap summary =====
top_final = set(top_final_df["guide_key"])
top_on_target = set(top_on_target_df["guide_key"])
top_off_target = set(top_off_target_df["guide_key"])
top_conservation = set(top_conservation_df["guide_key"])
top_consensus = set(top_consensus_df["guide_key"])

overlap_summary = pd.DataFrame({
    "comparison": [
        "Final vs On-target",
        "Final vs Off-target",
        "Final vs Conservation",
        "Final vs Consensus",
    ],
    "top10_overlap_count": [
        len(top_final & top_on_target),
        len(top_final & top_off_target),
        len(top_final & top_conservation),
        len(top_final & top_consensus),
    ]
})

overlap_output = RESULTS_DIR / "cross_model_top10_overlap_summary.csv"
overlap_summary.to_csv(overlap_output, index=False)

# ===== Print saved files =====
print("\nSaved files:")
print("-", full_ranked_output)
print("-", RESULTS_DIR / "top10_on_target_model.csv")
print("-", RESULTS_DIR / "top10_off_target_model.csv")
print("-", RESULTS_DIR / "top10_conservation_model.csv")
print("-", RESULTS_DIR / "top10_consensus_model.csv")
print("-", RESULTS_DIR / "top10_final_composite_model.csv")
print("-", comparison_output)
print("-", dropped_output)
print("-", overlap_output)
# ===== Figure: Cross-model agreement analysis =====
import matplotlib.pyplot as plt

figure_df = overlap_summary.copy()

label_map = {
    "Final vs On-target": "On-target-only\n(Activity baseline)",
    "Final vs Off-target": "Off-target-only\n(Specificity baseline)",
    "Final vs Conservation": "Conservation-only\n(Strain-awareness proxy)",
    "Final vs Consensus": "Consensus (Avg.)\n(Ensemble baseline)",
}

figure_df["plot_label"] = figure_df["comparison"].map(label_map)

plt.figure(figsize=(11, 7))

bars = plt.bar(
    figure_df["plot_label"],
    figure_df["top10_overlap_count"]
)

plt.title("Cross-model Agreement Analysis", fontsize=20, fontweight="bold", pad=18)
plt.suptitle("Top 10 Guide Overlap with the Proposed Composite Model (40/30/30)", fontsize=13, y=0.93)

plt.xlabel("Model Comparison", fontsize=13)
plt.ylabel("Overlap Count (Top 10)", fontsize=13)
plt.ylim(0, 11)
plt.yticks(range(0, 11, 2))
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.gca().set_axisbelow(True)

for bar, value in zip(bars, figure_df["top10_overlap_count"]):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.08,
        str(value),
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold"
    )

plt.tight_layout(rect=[0, 0, 1, 0.90])

figure_path = RESULTS_DIR / "Figure_CrossModel_Agreement_Analysis.png"
plt.savefig(figure_path, dpi=300, bbox_inches="tight")
plt.close()

print("-", figure_path)