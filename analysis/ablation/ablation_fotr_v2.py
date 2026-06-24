"""
Ablation analysis for FOTR v2.

Purpose:
- Test how much each scoring component affects guide ranking.
- Remove one component at a time:
    activity
    specificity
    conservation
    RNA accessibility
    functional context
    risk penalties
- Compare each ablated ranking against full FOTR v2 ranking.

Input:
    results_fotr_v2/all_guides_fotr_v2_functional_context_recommended.csv

Outputs:
    results_ablation/fotr_v2_ablation_all_scores.csv
    results_ablation/fotr_v2_ablation_summary.csv
    results_ablation/fotr_v2_ablation_top10_overlap_by_gene.csv
    results_ablation/fotr_v2_ablation_rank_shift_examples.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np


INPUT_PATH = Path("results_fotr_v2/all_guides_fotr_v2_functional_context_recommended.csv")
OUTPUT_DIR = Path("results_ablation")

OUTPUT_ALL = OUTPUT_DIR / "fotr_v2_ablation_all_scores.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "fotr_v2_ablation_summary.csv"
OUTPUT_TOP10 = OUTPUT_DIR / "fotr_v2_ablation_top10_overlap_by_gene.csv"
OUTPUT_EXAMPLES = OUTPUT_DIR / "fotr_v2_ablation_rank_shift_examples.csv"


COMPONENTS = {
    "activity": "activity_component",
    "specificity": "specificity_component",
    "conservation": "conservation_component",
    "rna_accessibility": "rna_accessibility_component",
    "functional_context": "target_functional_context_component_v2",
}


FULL_WEIGHTS = {
    "activity": 0.30,
    "specificity": 0.20,
    "conservation": 0.20,
    "rna_accessibility": 0.15,
    "functional_context": 0.15,
}


def safe_numeric(df, col, default=0.0):
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(default)
    return df


def compute_score(df, weights, include_risk_penalty=True):
    score = np.zeros(len(df))

    for component_name, weight in weights.items():
        col = COMPONENTS[component_name]
        score += weight * df[col]

    score = 100 * score

    if include_risk_penalty:
        score -= 10 * df["binding_risk_scaled_v2"]
        score -= 10 * df["rna_structure_risk_scaled_v2"]

    return score.round(4)


def top_n_set(group, rank_col, n=10):
    return set(group.sort_values(rank_col).head(n)["guide_key"])


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    required = {
        "gene",
        "position",
        "spacer",
        "pam",
        "fotr_v2_priority_score",
        "activity_component",
        "specificity_component",
        "conservation_component",
        "rna_accessibility_component",
        "target_functional_context_component_v2",
        "binding_risk_scaled_v2",
        "rna_structure_risk_scaled_v2",
        "rna_accessibility_class",
        "fotr_v2_recommendation_status",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    numeric_cols = [
        "fotr_v2_priority_score",
        "activity_component",
        "specificity_component",
        "conservation_component",
        "rna_accessibility_component",
        "target_functional_context_component_v2",
        "binding_risk_scaled_v2",
        "rna_structure_risk_scaled_v2",
    ]

    for col in numeric_cols:
        df = safe_numeric(df, col)

    df["guide_key"] = (
        df["gene"].astype(str)
        + "|"
        + df["position"].astype(str)
        + "|"
        + df["spacer"].astype(str)
        + "|"
        + df["pam"].astype(str)
    )

    # Full score and rank
    df["full_fotr_v2_score"] = df["fotr_v2_priority_score"]
    df["full_rank_within_gene"] = (
        df.groupby("gene")["full_fotr_v2_score"]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    ablation_names = []

    # Remove each positive component, then renormalize remaining positive weights to sum to 1.
    for removed_component in COMPONENTS:
        remaining = {k: v for k, v in FULL_WEIGHTS.items() if k != removed_component}
        remaining_sum = sum(remaining.values())
        renorm_weights = {k: v / remaining_sum for k, v in remaining.items()}

        score_col = f"score_without_{removed_component}"
        rank_col = f"rank_without_{removed_component}"
        delta_col = f"rank_delta_without_{removed_component}"

        df[score_col] = compute_score(df, renorm_weights, include_risk_penalty=True)
        df[rank_col] = (
            df.groupby("gene")[score_col]
            .rank(method="min", ascending=False)
            .astype(int)
        )
        df[delta_col] = df["full_rank_within_gene"] - df[rank_col]
        ablation_names.append(removed_component)

    # Remove risk penalty only
    score_col = "score_without_risk_penalties"
    rank_col = "rank_without_risk_penalties"
    delta_col = "rank_delta_without_risk_penalties"

    df[score_col] = compute_score(df, FULL_WEIGHTS, include_risk_penalty=False)
    df[rank_col] = (
        df.groupby("gene")[score_col]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    df[delta_col] = df["full_rank_within_gene"] - df[rank_col]
    ablation_names.append("risk_penalties")

    df.to_csv(OUTPUT_ALL, index=False)

    # Summary metrics
    summary_rows = []
    for name in ablation_names:
        score_col = f"score_without_{name}"
        rank_col = f"rank_without_{name}"
        delta_col = f"rank_delta_without_{name}"

        for gene, g in df.groupby("gene"):
            full_top = set(g.sort_values("full_rank_within_gene").head(10)["guide_key"])
            ablated_top = set(g.sort_values(rank_col).head(10)["guide_key"])
            overlap = len(full_top & ablated_top)
            jaccard = overlap / len(full_top | ablated_top) if len(full_top | ablated_top) else 0

            summary_rows.append({
                "gene": gene,
                "ablation": f"without_{name}",
                "mean_absolute_rank_shift": g[delta_col].abs().mean(),
                "max_absolute_rank_shift": g[delta_col].abs().max(),
                "top10_overlap_count": overlap,
                "top10_jaccard": round(jaccard, 4),
                "spearman_rank_correlation": round(
                    g["full_rank_within_gene"].corr(g[rank_col], method="spearman"),
                    4,
                ),
            })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_SUMMARY, index=False)

    top10 = summary[
        ["gene", "ablation", "top10_overlap_count", "top10_jaccard", "spearman_rank_correlation"]
    ].copy()
    top10.to_csv(OUTPUT_TOP10, index=False)

    # Most shifted examples across all ablations
    example_rows = []
    base_cols = [
        "gene",
        "position",
        "spacer",
        "pam",
        "rna_accessibility_class",
        "fotr_v2_recommendation_status",
        "full_fotr_v2_score",
        "full_rank_within_gene",
    ]

    for name in ablation_names:
        rank_col = f"rank_without_{name}"
        delta_col = f"rank_delta_without_{name}"
        temp = df[base_cols + [rank_col, delta_col]].copy()
        temp["ablation"] = f"without_{name}"
        temp["absolute_rank_shift"] = temp[delta_col].abs()
        temp = temp.sort_values("absolute_rank_shift", ascending=False).head(10)
        example_rows.append(temp)

    examples = pd.concat(example_rows, ignore_index=True)
    examples.to_csv(OUTPUT_EXAMPLES, index=False)

    print(f"Wrote: {OUTPUT_ALL}")
    print(f"Wrote: {OUTPUT_SUMMARY}")
    print(f"Wrote: {OUTPUT_TOP10}")
    print(f"Wrote: {OUTPUT_EXAMPLES}")

    print("\nAblation summary:")
    print(summary.sort_values(["gene", "mean_absolute_rank_shift"], ascending=[True, False]).to_string(index=False))

    print("\nLargest rank-shift examples:")
    print(examples.head(20).to_string(index=False))


if __name__ == "__main__":
    main()