from pathlib import Path
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]

TABLE_DIR = ROOT / "results_journal_extension/tables"
FIG_DIR = ROOT / "results_journal_extension/figures_bw"
REPORT_DIR = ROOT / "results_journal_extension/reports"

FIG_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

ALL20_FILE = TABLE_DIR / "journal_table_all20_integrated_fotr_guides.csv"
TOP2_FILE = TABLE_DIR / "journal_table_top2_guides_per_gene.csv"
STATUS_FILE = TABLE_DIR / "journal_table_requirement_completion_status.csv"

if not ALL20_FILE.exists():
    raise FileNotFoundError(
        f"Missing {ALL20_FILE}. Run this first:\n"
        f"python analysis/journal_extension/create_journal_ready_tables.py"
    )

all20 = pd.read_csv(ALL20_FILE)

if TOP2_FILE.exists():
    top2 = pd.read_csv(TOP2_FILE)
else:
    top2 = all20[all20["final_rank_within_gene"] <= 2].copy()

status = pd.read_csv(STATUS_FILE) if STATUS_FILE.exists() else None

gene_order = ["blaKPC", "blaNDM1", "mcr1", "mecA"]

all20["gene_order"] = all20["target_gene"].apply(
    lambda x: gene_order.index(x) if x in gene_order else 999
)
top2["gene_order"] = top2["target_gene"].apply(
    lambda x: gene_order.index(x) if x in gene_order else 999
)

def extract_original_rank(guide_id):
    m = re.search(r"top(\d+)", str(guide_id))
    return int(m.group(1)) if m else np.nan

all20["original_rank"] = all20["guide_id"].apply(extract_original_rank)

# -----------------------------
# Figure 1: Top two guides by conserved burden
# -----------------------------
fig1_df = top2.sort_values(["gene_order", "final_rank_within_gene"]).copy()
fig1_df["label"] = fig1_df["target_gene"] + " | " + fig1_df["guide_id"]

fig, ax = plt.subplots(figsize=(12, 8))
ypos = np.arange(len(fig1_df))

bars = ax.barh(
    ypos,
    fig1_df["conserved_offtarget_burden_score"],
    color="white",
    edgecolor="black",
    linewidth=1.2,
)

hatches = {
    "blaKPC": "///",
    "blaNDM1": "xxx",
    "mcr1": "...",
    "mecA": "\\\\\\",
}

for bar, gene in zip(bars, fig1_df["target_gene"]):
    bar.set_hatch(hatches.get(gene, ""))

max_x = max(float(fig1_df["conserved_offtarget_burden_score"].max()), 1.0)

for i, row in enumerate(fig1_df.itertuples(index=False)):
    ax.text(
        getattr(row, "conserved_offtarget_burden_score") + max_x * 0.03,
        i,
        f"hits={int(getattr(row, 'total_offtarget_hits'))}, acc={int(getattr(row, 'accessions_with_any_hit'))}",
        va="center",
        fontsize=8,
    )

ax.set_yticks(ypos)
ax.set_yticklabels(fig1_df["label"])
ax.invert_yaxis()
ax.set_xlabel("Conserved off-target burden score")
ax.set_title("Figure 1. Final top guides after integrated FOTR screening")
ax.grid(axis="x", linestyle="--", linewidth=0.5, color="black", alpha=0.4)
ax.set_xlim(0, max_x * 1.4)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure1_top_guides_conserved_burden_bw.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "figure1_top_guides_conserved_burden_bw.pdf", bbox_inches="tight")
plt.close()

# -----------------------------
# Figure 2: Original rank vs final integrated rank
# -----------------------------
rank_df = all20.sort_values(["gene_order", "original_rank"]).reset_index(drop=True)
rank_df["xpos"] = range(len(rank_df))

fig, ax = plt.subplots(figsize=(15, 8))

for _, row in rank_df.iterrows():
    x = row["xpos"]
    y1 = row["original_rank"]
    y2 = row["final_rank_within_gene"]

    ax.plot([x, x], [y1, y2], color="black", linewidth=1.1)
    ax.scatter(x, y1, s=60, marker="o", facecolors="white", edgecolors="black", zorder=3)
    ax.scatter(x, y2, s=45, marker="s", facecolors="black", edgecolors="black", zorder=3)

for gene in gene_order:
    subset = rank_df[rank_df["target_gene"] == gene]
    if len(subset):
        xmin = subset["xpos"].min()
        xmax = subset["xpos"].max()
        ax.text((xmin + xmax) / 2, 0.55, gene, ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.axvline(xmax + 0.5, color="black", linestyle=":", linewidth=0.8)

ax.set_xticks(rank_df["xpos"])
ax.set_xticklabels(rank_df["guide_id"], rotation=90, fontsize=8)
ax.set_ylabel("Rank within gene, 1 = best")
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
plt.savefig(FIG_DIR / "figure2_rank_shift_original_vs_final_bw.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "figure2_rank_shift_original_vs_final_bw.pdf", bbox_inches="tight")
plt.close()

# -----------------------------
# Figure 3: All-20 burden heatmap
# -----------------------------
heat_cols = [
    "total_offtarget_hits",
    "accessions_with_any_hit",
    "functional_burden_score",
    "mobile_context_burden_score",
    "conserved_offtarget_burden_score",
]

heat_df = all20.sort_values(["gene_order", "final_rank_within_gene"]).copy()
matrix = heat_df[heat_cols].astype(float)

norm_matrix = matrix.copy()
for col in heat_cols:
    colmax = norm_matrix[col].max()
    colmin = norm_matrix[col].min()
    if colmax == colmin:
        norm_matrix[col] = 0
    else:
        norm_matrix[col] = (norm_matrix[col] - colmin) / (colmax - colmin)

fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(norm_matrix.values, aspect="auto", cmap="Greys")

ax.set_xticks(np.arange(len(heat_cols)))
ax.set_xticklabels(
    ["Hits", "Accessions", "Functional\nburden", "Mobile\nburden", "Conserved\nburden"]
)

ax.set_yticks(np.arange(len(heat_df)))
ax.set_yticklabels(heat_df["target_gene"] + " | " + heat_df["guide_id"], fontsize=7)

for i in range(norm_matrix.shape[0]):
    for j in range(norm_matrix.shape[1]):
        raw_val = matrix.iloc[i, j]
        ax.text(j, i, f"{raw_val:.0f}", ha="center", va="center", fontsize=6)

ax.set_title("Figure 3. Integrated FOTR burden profile across all 20 guides")
plt.colorbar(im, ax=ax, label="Column-normalized burden intensity")
plt.tight_layout()
plt.savefig(FIG_DIR / "figure3_all20_integrated_burden_heatmap_bw.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "figure3_all20_integrated_burden_heatmap_bw.pdf", bbox_inches="tight")
plt.close()

# -----------------------------
# Figure 4: Total off-target hits by target gene
# -----------------------------
gene_hits = all20.groupby("target_gene", as_index=False).agg(
    total_hits=("total_offtarget_hits", "sum"),
    total_accession_hits=("accessions_with_any_hit", "sum"),
)

gene_hits["gene_order"] = gene_hits["target_gene"].apply(
    lambda x: gene_order.index(x) if x in gene_order else 999
)
gene_hits = gene_hits.sort_values("gene_order")

fig, ax = plt.subplots(figsize=(8, 6))
bars = ax.bar(
    gene_hits["target_gene"],
    gene_hits["total_hits"],
    color="white",
    edgecolor="black",
    linewidth=1.2,
)

for bar in bars:
    bar.set_hatch("///")

max_hits = max(float(gene_hits["total_hits"].max()), 1.0)

for i, row in enumerate(gene_hits.itertuples(index=False)):
    ax.text(
        i,
        getattr(row, "total_hits") + max_hits * 0.03,
        str(int(getattr(row, "total_hits"))),
        ha="center",
        va="bottom",
        fontsize=9,
    )

ax.set_ylabel("Total expanded off-target hit rows")
ax.set_xlabel("Target gene")
ax.set_title("Figure 4. Expanded off-target hit burden by AMR target")
ax.grid(axis="y", linestyle="--", linewidth=0.5, color="black", alpha=0.4)
ax.set_ylim(0, max_hits * 1.2)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure4_total_hits_by_target_gene_bw.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "figure4_total_hits_by_target_gene_bw.pdf", bbox_inches="tight")
plt.close()

# -----------------------------
# Figure 5: Completion status
# -----------------------------
if status is not None:
    plot_status = status.copy()

    def score_status(x):
        x = str(x)
        if "Complete computationally" in x:
            return 1.0
        if "Planning complete" in x:
            return 0.7
        return 0.0

    plot_status["completion_score"] = plot_status["journal_status"].apply(score_status)

    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(plot_status))
    bars = ax.barh(
        y,
        plot_status["completion_score"],
        color="white",
        edgecolor="black",
        linewidth=1.2,
    )

    for bar, stat in zip(bars, plot_status["journal_status"]):
        if "Planning" in str(stat):
            bar.set_hatch("...")
        elif "Not done" in str(stat):
            bar.set_hatch("xxx")
        else:
            bar.set_hatch("///")

    ax.set_yticks(y)
    ax.set_yticklabels(plot_status["analysis_layer"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Completion score")
    ax.set_title("Figure 5. Computational completion status for journal extension")
    ax.grid(axis="x", linestyle="--", linewidth=0.5, color="black", alpha=0.4)

    for i, row in plot_status.iterrows():
        ax.text(row["completion_score"] + 0.02, i, row["journal_status"], va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "figure5_journal_completion_status_bw.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIG_DIR / "figure5_journal_completion_status_bw.pdf", bbox_inches="tight")
    plt.close()

summary = []
summary.append("Black-and-white journal figures created.")
summary.append("")
summary.append("Figures:")
for f in sorted(FIG_DIR.glob("*.png")):
    summary.append(f"- {f}")
summary.append("")
summary.append("Recommended use:")
summary.append("- BIBM main track: use Figures 1, 2, and possibly 4.")
summary.append("- Journal extension: use Figures 1 to 5.")
summary.append("- Supplement: include full tables and extra heatmaps.")

(REPORT_DIR / "journal_figure_generation_summary.txt").write_text("\n".join(summary))

print("Journal black-and-white figures created.")
for f in sorted(FIG_DIR.glob("*.png")):
    print("Wrote:", f)
print("Wrote:", REPORT_DIR / "journal_figure_generation_summary.txt")
