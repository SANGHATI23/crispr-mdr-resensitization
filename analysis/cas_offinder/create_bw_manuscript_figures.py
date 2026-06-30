from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]

input_file = ROOT / "results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/final_requirement4_complete/requirement4_final_integrated_conserved_offtarget_table_ALL20.csv"
out_dir = ROOT / "results_cas_offinder/expanded_panel/manuscript_figures_bw"
out_dir.mkdir(parents=True, exist_ok=True)

if not input_file.exists():
    raise FileNotFoundError(f"Input file not found: {input_file}")

df = pd.read_csv(input_file)

gene_col = "target_gene"
guide_col = "guide_id"
hits_col = "total_offtarget_hit_rows"
acc_col = "accessions_with_any_hit"
func_col = "total_db_refined_burden"
mobile_col = "requirement3_mobile_context_burden"
cons_col = "conserved_functional_burden_score"
decision_col = "final_req4_integrated_decision"
final_rank_col = "final_req4_rank_within_gene"

required = [
    gene_col,
    guide_col,
    hits_col,
    acc_col,
    func_col,
    mobile_col,
    cons_col,
    decision_col,
    final_rank_col,
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}\nAvailable columns: {list(df.columns)}")

df = df.copy()

def extract_original_rank(guide_id):
    m = re.search(r"top(\d+)", str(guide_id))
    return int(m.group(1)) if m else None

df["original_rank"] = df[guide_col].apply(extract_original_rank)

gene_order = ["blaKPC", "blaNDM1", "mcr1", "mecA"]
df["gene_order"] = df[gene_col].apply(lambda x: gene_order.index(x) if x in gene_order else 999)

df["final_rank"] = df[final_rank_col]

# -----------------------------
# Figure 1: Top 2 guides per gene
# -----------------------------
top2 = df[df["final_rank"] <= 2].copy()
top2 = top2.sort_values(["gene_order", "final_rank"]).reset_index(drop=True)
top2["label"] = top2[gene_col] + " | " + top2[guide_col]

hatch_map = {
    "blaKPC": "///",
    "blaNDM1": "xxx",
    "mcr1": "...",
    "mecA": "\\\\\\"
}

fig, ax = plt.subplots(figsize=(12, 8))

ypos = list(range(len(top2)))
bars = ax.barh(
    ypos,
    top2[cons_col],
    color="white",
    edgecolor="black",
    linewidth=1.2
)

for bar, gene in zip(bars, top2[gene_col]):
    bar.set_hatch(hatch_map.get(gene, ""))

max_val = max(top2[cons_col].max(), 1)

for i, row in top2.iterrows():
    txt = f"hits={int(row[hits_col])}, acc={int(row[acc_col])}"
    ax.text(
        row[cons_col] + max_val * 0.03,
        i,
        txt,
        va="center",
        fontsize=8
    )

ax.set_yticks(ypos)
ax.set_yticklabels(top2["label"])
ax.invert_yaxis()
ax.set_xlabel("Conserved off-target burden score (lower is better)")
ax.set_title("Figure 1. Final top two guides per gene after integrated FOTR screening")
ax.grid(axis="x", linestyle="--", linewidth=0.5, color="black", alpha=0.4)

ax.set_xlim(0, max_val * 1.35)

plt.tight_layout()
fig1_png = out_dir / "figure1_top2_guides_conserved_burden_bw.png"
fig1_pdf = out_dir / "figure1_top2_guides_conserved_burden_bw.pdf"
plt.savefig(fig1_png, dpi=300, bbox_inches="tight")
plt.savefig(fig1_pdf, bbox_inches="tight")
plt.close()

# -----------------------------
# Figure 2: Original rank vs final integrated rank
# -----------------------------
rank_df = df.sort_values(["gene_order", "original_rank"]).reset_index(drop=True)
rank_df["xpos"] = range(len(rank_df))

fig, ax = plt.subplots(figsize=(15, 8))

for _, row in rank_df.iterrows():
    x = row["xpos"]
    y1 = row["original_rank"]
    y2 = row["final_rank"]

    ax.plot([x, x], [y1, y2], color="black", linewidth=1.2)
    ax.scatter(x, y1, s=60, marker="o", facecolors="white", edgecolors="black", zorder=3)
    ax.scatter(x, y2, s=45, marker="s", facecolors="black", edgecolors="black", zorder=3)

for gene in gene_order:
    subset = rank_df[rank_df[gene_col] == gene]
    if not subset.empty:
        xmin = subset["xpos"].min()
        xmax = subset["xpos"].max()
        xmid = (xmin + xmax) / 2
        ax.text(xmid, 0.55, gene, ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.axvline(xmax + 0.5, color="black", linestyle=":", linewidth=0.8)

ax.set_xticks(rank_df["xpos"])
ax.set_xticklabels(rank_df[guide_col], rotation=90, fontsize=8)
ax.set_ylabel("Rank within gene (1 = best)")
ax.set_title("Figure 2. Original risk-aware rank versus final integrated rank")
ax.invert_yaxis()
ax.set_ylim(5.5, 0.4)
ax.grid(axis="y", linestyle="--", linewidth=0.5, color="black", alpha=0.4)

legend_elements = [
    Line2D([0], [0], marker="o", color="black", markerfacecolor="white",
           markersize=7, linestyle="None", label="Original risk-aware rank"),
    Line2D([0], [0], marker="s", color="black", markerfacecolor="black",
           markersize=6, linestyle="None", label="Final integrated rank"),
]
ax.legend(handles=legend_elements, loc="upper right", frameon=True)

plt.tight_layout()
fig2_png = out_dir / "figure2_rank_shift_original_vs_final_bw.png"
fig2_pdf = out_dir / "figure2_rank_shift_original_vs_final_bw.pdf"
plt.savefig(fig2_png, dpi=300, bbox_inches="tight")
plt.savefig(fig2_pdf, bbox_inches="tight")
plt.close()

summary_file = out_dir / "figure_generation_summary.txt"
summary_text = f"""
Black-and-white manuscript figures created successfully.

Input:
{input_file}

Outputs:
{fig1_png}
{fig1_pdf}
{fig2_png}
{fig2_pdf}

Figure 1:
Top two guides per gene after integrated FOTR screening, plotted by conserved off-target burden.

Figure 2:
Original risk-aware rank versus final integrated rank within each gene.

Interpretation:
- Lower conserved off-target burden is better.
- Rank shifts show the effect of adding genome-wide off-target hits, functional annotation, mobile-context mapping, and strain-level off-target conservation.
- blaKPC_riskaware_top5 and blaNDM1_riskaware_top1 are expected to appear as the cleanest computational candidates.
"""
summary_file.write_text(summary_text.strip())

print("Black-and-white manuscript figures created successfully.")
print("Wrote:", fig1_png)
print("Wrote:", fig1_pdf)
print("Wrote:", fig2_png)
print("Wrote:", fig2_pdf)
print("Wrote:", summary_file)
