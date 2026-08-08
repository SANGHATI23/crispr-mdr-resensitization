"""
FOTR v2 target-functional-context scoring for CRISPR-MDR.

Purpose:
- Upgrade FOTR v1 by replacing neutral placeholder functional/context terms
  with target-side functional annotation v1.
- This is still not full genome-wide off-target functional annotation.
- It is a transparent FOTR v2 intermediate step for BIBM development.

Input:
    results_functional_annotation/all_guides_functional_context_v1.csv

Outputs:
    results_fotr_v2/all_guides_fotr_v2_functional_context.csv
    results_fotr_v2/top5_fotr_v2_by_gene.csv
    results_fotr_v2/fotr_v1_vs_v2_rank_change.csv
    results_fotr_v2/fotr_v2_summary_by_gene.csv
"""

from pathlib import Path
import pandas as pd
import numpy as np


INPUT_PATH = Path("results_functional_annotation/all_guides_functional_context_v1.csv")
OUTPUT_DIR = Path("results_fotr_v2")

OUTPUT_ALL = OUTPUT_DIR / "all_guides_fotr_v2_functional_context.csv"
OUTPUT_TOP5 = OUTPUT_DIR / "top5_fotr_v2_by_gene.csv"
OUTPUT_RANK_CHANGE = OUTPUT_DIR / "fotr_v1_vs_v2_rank_change.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "fotr_v2_summary_by_gene.csv"


def safe_numeric(series, default=0.0):
    return pd.to_numeric(series, errors="coerce").fillna(default)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    required_cols = {
        "gene",
        "spacer",
        "pam",
        "position",
        "activity_component",
        "specificity_component",
        "conservation_component",
        "rna_accessibility_component",
        "binding_risk",
        "rna_structure_risk",
        "functional_severity_v1",
        "target_context_weight_v1",
        "fotr_priority_score",
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Normalize numeric inputs
    df["activity_component"] = safe_numeric(df["activity_component"])
    df["specificity_component"] = safe_numeric(df["specificity_component"])
    df["conservation_component"] = safe_numeric(df["conservation_component"])
    df["rna_accessibility_component"] = safe_numeric(df["rna_accessibility_component"])
    df["binding_risk"] = safe_numeric(df["binding_risk"])
    df["rna_structure_risk"] = safe_numeric(df["rna_structure_risk"])
    df["functional_severity_v1"] = safe_numeric(df["functional_severity_v1"], default=0.75)
    df["target_context_weight_v1"] = safe_numeric(df["target_context_weight_v1"], default=0.75)

    # Convert risk values to bounded penalties.
    # binding_risk in current table often equals 0 for perfect specificity rows.
    # rna_structure_risk can be larger, so scale to 0-1.
    df["binding_risk_scaled_v2"] = np.clip(df["binding_risk"] / 100.0, 0, 1)
    df["rna_structure_risk_scaled_v2"] = np.clip(df["rna_structure_risk"] / 100.0, 0, 1)

    # Target functional benefit: high if the guide hits functionally important target context.
    df["target_functional_context_component_v2"] = (
        df["functional_severity_v1"] * df["target_context_weight_v1"]
    ).clip(0, 1)

    # FOTR v2 priority score:
    # Positive components reward activity, specificity, conservation, RNA accessibility,
    # and target functional/context importance.
    # Negative components penalize binding risk and RNA-structure risk.
    #
    # Weights are intentionally transparent and will later be tested through ablation/sensitivity.
    df["fotr_v2_priority_score"] = (
        100
        * (
            0.30 * df["activity_component"]
            + 0.20 * df["specificity_component"]
            + 0.20 * df["conservation_component"]
            + 0.15 * df["rna_accessibility_component"]
            + 0.15 * df["target_functional_context_component_v2"]
        )
        - 10 * df["binding_risk_scaled_v2"]
        - 10 * df["rna_structure_risk_scaled_v2"]
    )

    df["fotr_v2_priority_score"] = df["fotr_v2_priority_score"].round(4)

    # Keep old FOTR v1 score as comparison
    df["fotr_v1_priority_score"] = safe_numeric(df["fotr_priority_score"])

    # Rank within gene for v1 and v2
    df["fotr_v1_rank_within_gene"] = (
        df.groupby("gene")["fotr_v1_priority_score"]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    df["fotr_v2_rank_within_gene"] = (
        df.groupby("gene")["fotr_v2_priority_score"]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    df["rank_change_v1_minus_v2"] = (
        df["fotr_v1_rank_within_gene"] - df["fotr_v2_rank_within_gene"]
    )

    df["rank_change_direction"] = np.where(
        df["rank_change_v1_minus_v2"] > 0,
        "improved_in_v2",
        np.where(
            df["rank_change_v1_minus_v2"] < 0,
            "dropped_in_v2",
            "unchanged",
        ),
    )

    # Save full table
    df.to_csv(OUTPUT_ALL, index=False)

    # Top 5 by gene
    top5 = (
        df.sort_values(["gene", "fotr_v2_priority_score"], ascending=[True, False])
        .groupby("gene")
        .head(5)
    )

    top5_cols = [
        "gene",
        "position",
        "spacer",
        "pam",
        "rna_accessibility_class",
        "target_region_bin",
        "functional_severity_v1",
        "target_context_weight_v1",
        "target_functional_context_component_v2",
        "fotr_v1_priority_score",
        "fotr_v2_priority_score",
        "fotr_v1_rank_within_gene",
        "fotr_v2_rank_within_gene",
        "rank_change_v1_minus_v2",
        "rank_change_direction",
        "crisot_available",
        "crisot_wt_specificity",
    ]

    top5[top5_cols].to_csv(OUTPUT_TOP5, index=False)

    # Rank change table
    rank_change_cols = [
        "gene",
        "position",
        "spacer",
        "pam",
        "rna_accessibility_class",
        "target_region_bin",
        "fotr_v1_priority_score",
        "fotr_v2_priority_score",
        "fotr_v1_rank_within_gene",
        "fotr_v2_rank_within_gene",
        "rank_change_v1_minus_v2",
        "rank_change_direction",
    ]

    rank_change = df[rank_change_cols].sort_values(
        ["gene", "rank_change_v1_minus_v2"],
        ascending=[True, False],
    )
    rank_change.to_csv(OUTPUT_RANK_CHANGE, index=False)

    # Summary
    summary = (
        df.groupby("gene")
        .agg(
            guide_rows=("spacer", "count"),
            unique_spacers=("spacer", "nunique"),
            mean_fotr_v1=("fotr_v1_priority_score", "mean"),
            mean_fotr_v2=("fotr_v2_priority_score", "mean"),
            max_fotr_v2=("fotr_v2_priority_score", "max"),
            accessible_guides=("rna_accessibility_class", lambda x: (x == "Accessible").sum()),
            early_coding_guides=("target_region_bin", lambda x: (x == "early_coding_region").sum()),
            mean_functional_context_component_v2=(
                "target_functional_context_component_v2",
                "mean",
            ),
        )
        .reset_index()
    )

    summary.to_csv(OUTPUT_SUMMARY, index=False)

    print(f"Wrote: {OUTPUT_ALL}")
    print(f"Wrote: {OUTPUT_TOP5}")
    print(f"Wrote: {OUTPUT_RANK_CHANGE}")
    print(f"Wrote: {OUTPUT_SUMMARY}")

    print("\nTop 5 FOTR v2 guides by gene:")
    print(top5[top5_cols].to_string(index=False))

    print("\nSummary by gene:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()