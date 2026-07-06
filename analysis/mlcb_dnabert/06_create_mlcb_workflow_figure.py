
#!/usr/bin/env python3
"""
Create Figure 1 for the MLCB manuscript.

Figure 1 shows the full context-aware FOTR-CRISPR + DNABERT-2 workflow:
AMR references -> guide generation -> 100-bp local context extraction ->
frozen DNABERT-2 embedding -> FOTR annotation features ->
weak functional severity label -> leakage-aware ablation/risk model ->
risk-aware guide ranking.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = PROJECT_ROOT / "results_mlcb" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PNG_OUT = FIGURE_DIR / "figure1_mlcb_workflow.png"
PDF_OUT = FIGURE_DIR / "figure1_mlcb_workflow.pdf"


def add_box(ax, x, y, w, h, title, subtitle=None):
    """Add rounded workflow box."""
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.4,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(box)

    ax.text(
        x + w / 2,
        y + h * 0.62,
        title,
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        wrap=True,
    )

    if subtitle:
        ax.text(
            x + w / 2,
            y + h * 0.32,
            subtitle,
            ha="center",
            va="center",
            fontsize=8.5,
            wrap=True,
        )


def add_arrow(ax, x1, y1, x2, y2):
    """Add workflow arrow."""
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=15,
        linewidth=1.3,
        color="black",
        shrinkA=4,
        shrinkB=4,
    )
    ax.add_patch(arrow)


def main():
    fig, ax = plt.subplots(figsize=(8.5, 10.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title
    ax.text(
        0.5,
        0.965,
        "Figure 1. Context-aware FOTR-CRISPR + DNABERT-2 workflow",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )

    # Workflow boxes
    x = 0.17
    w = 0.66
    h = 0.075

    boxes = [
        (
            0.85,
            "AMR genome / multistrain references",
            "AMR target FASTA files and guide-level candidate table",
        ),
        (
            0.735,
            "Guide generation",
            "Candidate CRISPR-AMR spacers for blaNDM1, blaKPC, mcr1, and mecA",
        ),
        (
            0.62,
            "100-bp local genomic context extraction",
            "Target-centered windows; reverse-strand matches reverse-complemented",
        ),
        (
            0.505,
            "Frozen DNABERT-2 embedding",
            "100-bp context windows represented as 768-dimensional vectors",
        ),
        (
            0.39,
            "FOTR annotation features",
            "Guide-level biological and scoring-derived annotations",
        ),
        (
            0.275,
            "Weak functional severity label",
            "Annotation-derived computational label, not experimental ground truth",
        ),
        (
            0.16,
            "Leakage-aware ablation / risk model",
            "Score-, rank-, label-, and coordinate-derived features removed",
        ),
        (
            0.045,
            "Risk-aware guide ranking",
            "Candidate guides prioritized under functional off-target risk penalty",
        ),
    ]

    for y, title, subtitle in boxes:
        add_box(ax, x, y, w, h, title, subtitle)

    # Arrows between boxes
    for i in range(len(boxes) - 1):
        y_current = boxes[i][0]
        y_next = boxes[i + 1][0]
        add_arrow(ax, 0.5, y_current, 0.5, y_next + h)

    # Side note for ML contribution
    note_x = 0.04
    note_y = 0.47
    note_w = 0.11
    note_h = 0.18

    note = FancyBboxPatch(
        (note_x, note_y),
        note_w,
        note_h,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=1.0,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(note)

    ax.text(
        note_x + note_w / 2,
        note_y + note_h / 2,
        "ML layer:\nlocal-context\nfoundation-model\nrepresentations",
        ha="center",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        wrap=True,
    )

    add_arrow(ax, note_x + note_w, note_y + note_h / 2, x, 0.505 + h / 2)

    # Side note for leakage control
    note2_x = 0.85
    note2_y = 0.205
    note2_w = 0.11
    note2_h = 0.17

    note2 = FancyBboxPatch(
        (note2_x, note2_y),
        note2_w,
        note2_h,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=1.0,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(note2)

    ax.text(
        note2_x + note2_w / 2,
        note2_y + note2_h / 2,
        "Leakage control:\nremove direct\nweak-label\nconstruction\nfeatures",
        ha="center",
        va="center",
        fontsize=8.0,
        fontweight="bold",
        wrap=True,
    )

    add_arrow(ax, x + w, 0.16 + h / 2, note2_x, note2_y + note2_h / 2)

    # Footer
    ax.text(
        0.5,
        0.015,
        "Final manuscript analysis uses 100-bp local-context DNABERT-2 inputs and leakage-aware ablation.",
        ha="center",
        va="center",
        fontsize=8.5,
    )

    fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote: {PNG_OUT}")
    print(f"Wrote: {PDF_OUT}")


if __name__ == "__main__":
    main()
