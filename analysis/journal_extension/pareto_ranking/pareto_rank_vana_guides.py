#!/usr/bin/env python3

"""
pareto_rank_vana_guides.py

Purpose:
Create Pareto-front ranking for prospective unseen-gene vanA CRISPR guides.

This combines:
    1. First-pass unseen-gene priority score
    2. Escape-resilience score
    3. Simple on-target proxy
    4. GC acceptability
    5. Target-region priority

Inputs:
    results_journal_extension/unseen_genes/unseen_gene_candidate_guides_ranked_v1.csv
    results_journal_extension/escape_simulation/vanA_escape_summary_by_guide.csv

Outputs:
    results_journal_extension/pareto_ranking/vanA_pareto_ranked_guides.csv
    results_journal_extension/pareto_ranking/vanA_pareto_front_guides.csv
    results_journal_extension/pareto_ranking/vanA_pareto_summary.txt
"""

from pathlib import Path
import pandas as pd


GUIDE_INPUT = Path(
    "results_journal_extension/unseen_genes/"
    "unseen_gene_candidate_guides_ranked_v1.csv"
)

ESCAPE_INPUT = Path(
    "results_journal_extension/escape_simulation/"
    "vanA_escape_summary_by_guide.csv"
)

OUTPUT_DIR = Path("results_journal_extension/pareto_ranking")
OUTPUT_ALL = OUTPUT_DIR / "vanA_pareto_ranked_guides.csv"
OUTPUT_FRONT = OUTPUT_DIR / "vanA_pareto_front_guides.csv"
OUTPUT_TXT = OUTPUT_DIR / "vanA_pareto_summary.txt"


OBJECTIVE_COLUMNS = [
    "unseen_gene_priority_score_v1",
    "escape_resilience_score",
    "simple_on_target_proxy",
    "gc_objective_score",
    "target_region_objective_score",
]


def gc_objective(gc_class: str) -> float:
    """
    Convert GC class into a numeric objective.
    Higher is better.
    """
    gc_class = str(gc_class)

    if gc_class == "Acceptable_GC":
        return 1.0
    if gc_class in {"Moderate_GC", "Borderline_GC"}:
        return 0.7
    if gc_class in {"High_GC", "Low_GC"}:
        return 0.4
    if "Extreme" in gc_class:
        return 0.1
    return 0.5


def target_region_objective(region: str) -> float:
    """
    Convert target-region bin into a numeric objective.
    Higher is better.
    """
    region = str(region)

    if region == "early_coding_region":
        return 1.0
    if region == "middle_coding_region":
        return 0.7
    if region == "late_coding_region":
        return 0.4
    return 0.5


def dominates(row_a: pd.Series, row_b: pd.Series, objectives: list[str]) -> bool:
    """
    Return True if row_a Pareto-dominates row_b.

    A dominates B if:
        - A is at least as good as B in every objective.
        - A is strictly better than B in at least one objective.
    """
    at_least_equal_all = True
    strictly_better_any = False

    for col in objectives:
        a = row_a[col]
        b = row_b[col]

        if a < b:
            at_least_equal_all = False
            break

        if a > b:
            strictly_better_any = True

    return at_least_equal_all and strictly_better_any


def assign_pareto_front(df: pd.DataFrame, objectives: list[str]) -> pd.DataFrame:
    """
    Assign Pareto-front status.
    This first version identifies only the non-dominated front.
    """
    df = df.copy()
    is_pareto_front = []

    for i, row_i in df.iterrows():
        dominated = False

        for j, row_j in df.iterrows():
            if i == j:
                continue

            if dominates(row_j, row_i, objectives):
                dominated = True
                break

        is_pareto_front.append(not dominated)

    df["is_pareto_front"] = is_pareto_front
    return df


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not GUIDE_INPUT.exists():
        raise FileNotFoundError(f"Missing guide input: {GUIDE_INPUT}")

    if not ESCAPE_INPUT.exists():
        raise FileNotFoundError(f"Missing escape input: {ESCAPE_INPUT}")

    guides = pd.read_csv(GUIDE_INPUT)
    escape = pd.read_csv(ESCAPE_INPUT)

    if "gene" in guides.columns:
        guides = guides[guides["gene"].astype(str).str.lower() == "vana"].copy()

    if guides.empty:
        raise ValueError("No vanA guide rows found in guide input.")

    required_guide_cols = [
        "gene",
        "spacer",
        "pam",
        "strand",
        "position_1based",
        "gc_content",
        "gc_class",
        "target_region_bin",
        "simple_on_target_proxy",
        "unseen_gene_priority_score_v1",
        "unseen_gene_recommendation_v1",
    ]

    missing_guide_cols = [
        col for col in required_guide_cols if col not in guides.columns
    ]

    if missing_guide_cols:
        raise ValueError(f"Missing guide columns: {missing_guide_cols}")

    required_escape_cols = [
        "original_spacer",
        "pam",
        "escape_resilience_score",
        "escape_resilience_percentile",
        "escape_resilience_class",
    ]

    missing_escape_cols = [
        col for col in required_escape_cols if col not in escape.columns
    ]

    if missing_escape_cols:
        raise ValueError(f"Missing escape columns: {missing_escape_cols}")

    merged = guides.merge(
        escape[
            [
                "original_spacer",
                "pam",
                "escape_resilience_score",
                "escape_resilience_percentile",
                "escape_resilience_class",
                "escape_vulnerability_score",
                "total_single_mutants",
            ]
        ],
        left_on=["spacer", "pam"],
        right_on=["original_spacer", "pam"],
        how="left",
    )

    if merged["escape_resilience_score"].isna().any():
        missing_count = merged["escape_resilience_score"].isna().sum()
        raise ValueError(
            f"{missing_count} guides did not receive escape-resilience scores."
        )

    merged["gc_objective_score"] = merged["gc_class"].apply(gc_objective)
    merged["target_region_objective_score"] = merged["target_region_bin"].apply(
        target_region_objective
    )

    # Scale smaller objective values to a 0-100 style range for easier reading.
    merged["gc_objective_score_100"] = merged["gc_objective_score"] * 100
    merged["target_region_objective_score_100"] = (
        merged["target_region_objective_score"] * 100
    )

    pareto_df = assign_pareto_front(merged, OBJECTIVE_COLUMNS)

    # Helpful composite score for ordering only.
    # Pareto-front status remains the main multi-objective result.
    pareto_df["pareto_support_score"] = (
        0.35 * pareto_df["unseen_gene_priority_score_v1"]
        + 0.25 * pareto_df["escape_resilience_score"]
        + 0.20 * pareto_df["simple_on_target_proxy"]
        + 0.10 * pareto_df["gc_objective_score_100"]
        + 0.10 * pareto_df["target_region_objective_score_100"]
    ).round(2)

    pareto_df["pareto_recommendation"] = "Non_pareto_candidate"
    pareto_df.loc[
        pareto_df["is_pareto_front"], "pareto_recommendation"
    ] = "Pareto_front_candidate"

    pareto_df.loc[
        (pareto_df["is_pareto_front"])
        & (
            pareto_df["unseen_gene_recommendation_v1"]
            == "High_priority_unseen_candidate"
        )
        & (pareto_df["escape_resilience_class"] == "High_escape_resilience"),
        "pareto_recommendation",
    ] = "High_priority_pareto_candidate"

    pareto_df = pareto_df.sort_values(
        by=[
            "is_pareto_front",
            "pareto_support_score",
            "unseen_gene_priority_score_v1",
            "escape_resilience_score",
        ],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)

    pareto_df["pareto_rank"] = range(1, len(pareto_df) + 1)

    pareto_front = pareto_df[pareto_df["is_pareto_front"]].copy()

    output_columns = [
        "pareto_rank",
        "gene",
        "spacer",
        "pam",
        "strand",
        "position_1based",
        "gc_content",
        "gc_class",
        "target_region_bin",
        "simple_on_target_proxy",
        "unseen_gene_priority_score_v1",
        "unseen_gene_recommendation_v1",
        "escape_resilience_score",
        "escape_resilience_percentile",
        "escape_resilience_class",
        "gc_objective_score_100",
        "target_region_objective_score_100",
        "pareto_support_score",
        "is_pareto_front",
        "pareto_recommendation",
    ]

    pareto_df[output_columns].to_csv(OUTPUT_ALL, index=False)
    pareto_front[output_columns].to_csv(OUTPUT_FRONT, index=False)

    total_guides = len(pareto_df)
    pareto_count = len(pareto_front)
    high_priority_pareto_count = (
        pareto_df["pareto_recommendation"] == "High_priority_pareto_candidate"
    ).sum()

    top = pareto_df.iloc[0]

    with open(OUTPUT_TXT, "w") as f:
        f.write("vanA Pareto-Front Ranking Summary\n")
        f.write("=================================\n\n")
        f.write(f"Total vanA guides analyzed: {total_guides}\n")
        f.write(f"Pareto-front candidates: {pareto_count}\n")
        f.write(
            f"High-priority Pareto candidates: "
            f"{high_priority_pareto_count}\n\n"
        )

        f.write("Objectives optimized:\n")
        f.write("- unseen_gene_priority_score_v1: higher is better\n")
        f.write("- escape_resilience_score: higher is better\n")
        f.write("- simple_on_target_proxy: higher is better\n")
        f.write("- gc_objective_score: acceptable GC is better\n")
        f.write("- target_region_objective_score: early coding is better\n\n")

        f.write("Top Pareto-ranked guide:\n")
        f.write(f"- Spacer: {top['spacer']}\n")
        f.write(f"- PAM: {top['pam']}\n")
        f.write(f"- Position: {top['position_1based']}\n")
        f.write(f"- Target region: {top['target_region_bin']}\n")
        f.write(
            f"- Unseen-gene priority score v1: "
            f"{top['unseen_gene_priority_score_v1']}\n"
        )
        f.write(
            f"- Escape resilience score: "
            f"{top['escape_resilience_score']}\n"
        )
        f.write(
            f"- Escape resilience class: "
            f"{top['escape_resilience_class']}\n"
        )
        f.write(f"- Pareto support score: {top['pareto_support_score']}\n")
        f.write(f"- Pareto recommendation: {top['pareto_recommendation']}\n")

    print(f"Wrote: {OUTPUT_ALL}")
    print(f"Wrote: {OUTPUT_FRONT}")
    print(f"Wrote: {OUTPUT_TXT}")
    print()
    print("Pareto-front ranking completed.")
    print(f"Total vanA guides analyzed: {total_guides}")
    print(f"Pareto-front candidates: {pareto_count}")
    print(f"High-priority Pareto candidates: {high_priority_pareto_count}")


if __name__ == "__main__":
    main()