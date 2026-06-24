#!/usr/bin/env python3

"""
simulate_single_mutation_escape.py

Simulates first-pass single-nucleotide escape risk for vanA candidate guides.

Input:
    results_journal_extension/unseen_genes/unseen_gene_candidate_guides_ranked_v1.csv

Outputs:
    results_journal_extension/escape_simulation/vanA_escape_simulation_results.csv
    results_journal_extension/escape_simulation/vanA_escape_summary_by_guide.csv
    results_journal_extension/escape_simulation/vanA_escape_summary.txt
"""

from pathlib import Path
import pandas as pd


BASES = ["A", "C", "G", "T"]

INPUT_PATH = Path(
    "results_journal_extension/unseen_genes/unseen_gene_candidate_guides_ranked_v1.csv"
)

OUTPUT_DIR = Path("results_journal_extension/escape_simulation")
OUTPUT_MUTANTS = OUTPUT_DIR / "vanA_escape_simulation_results.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "vanA_escape_summary_by_guide.csv"
OUTPUT_TXT = OUTPUT_DIR / "vanA_escape_summary.txt"


def classify_escape_risk(position_1based: int) -> str:
    """
    First-pass escape-risk classification.

    Positions 1-7 are treated as highest-risk escape positions.
    Positions 8-12 are moderate-risk.
    Positions 13-20 are lower-risk.
    """
    if 1 <= position_1based <= 7:
        return "High_escape_risk"
    if 8 <= position_1based <= 12:
        return "Moderate_escape_risk"
    return "Lower_escape_risk"


def get_risk_weight(risk_class: str) -> float:
    if risk_class == "High_escape_risk":
        return 1.0
    if risk_class == "Moderate_escape_risk":
        return 0.5
    if risk_class == "Lower_escape_risk":
        return 0.1
    return 0.0


def generate_single_mutants(row: pd.Series) -> list[dict]:
    spacer = str(row["spacer"]).upper()
    gene = row.get("gene", "unknown_gene")

    mutant_rows = []

    for idx, original_base in enumerate(spacer):
        mutation_position = idx + 1

        for mutated_base in BASES:
            if mutated_base == original_base:
                continue

            mutant_spacer = spacer[:idx] + mutated_base + spacer[idx + 1:]
            escape_class = classify_escape_risk(mutation_position)

            mutant_rows.append(
                {
                    "gene": gene,
                    "original_spacer": spacer,
                    "pam": row.get("pam"),
                    "strand": row.get("strand"),
                    "original_position_1based": row.get("position_1based"),
                    "mutation_position_in_spacer": mutation_position,
                    "original_base": original_base,
                    "mutated_base": mutated_base,
                    "mutant_spacer": mutant_spacer,
                    "escape_risk_class": escape_class,
                    "escape_risk_weight": get_risk_weight(escape_class),
                    "unseen_gene_rank_v1": row.get("unseen_gene_rank_v1"),
                    "unseen_gene_priority_score_v1": row.get(
                        "unseen_gene_priority_score_v1"
                    ),
                    "unseen_gene_recommendation_v1": row.get(
                        "unseen_gene_recommendation_v1"
                    ),
                    "gc_content": row.get("gc_content"),
                    "gc_class": row.get("gc_class"),
                    "target_region_bin": row.get("target_region_bin"),
                    "simple_on_target_proxy": row.get("simple_on_target_proxy"),
                }
            )

    return mutant_rows


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    guides = pd.read_csv(INPUT_PATH)

    required_columns = ["spacer", "pam"]
    missing_columns = [col for col in required_columns if col not in guides.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if "gene" in guides.columns:
        guides = guides[guides["gene"].astype(str).str.lower() == "vana"].copy()

    all_mutants = []
    for _, row in guides.iterrows():
        all_mutants.extend(generate_single_mutants(row))

    mutants_df = pd.DataFrame(all_mutants)
    mutants_df.to_csv(OUTPUT_MUTANTS, index=False)

    summary = (
        mutants_df.groupby(["gene", "original_spacer", "pam"], dropna=False)
        .agg(
            total_single_mutants=("mutant_spacer", "count"),
            high_escape_mutants=(
                "escape_risk_class",
                lambda x: (x == "High_escape_risk").sum(),
            ),
            moderate_escape_mutants=(
                "escape_risk_class",
                lambda x: (x == "Moderate_escape_risk").sum(),
            ),
            low_escape_mutants=(
                "escape_risk_class",
                lambda x: (x == "Lower_escape_risk").sum(),
            ),
            escape_vulnerability_score=("escape_risk_weight", "mean"),
            unseen_gene_rank_v1=("unseen_gene_rank_v1", "first"),
            unseen_gene_priority_score_v1=(
                "unseen_gene_priority_score_v1",
                "first",
            ),
            unseen_gene_recommendation_v1=(
                "unseen_gene_recommendation_v1",
                "first",
            ),
            gc_content=("gc_content", "first"),
            gc_class=("gc_class", "first"),
            target_region_bin=("target_region_bin", "first"),
            simple_on_target_proxy=("simple_on_target_proxy", "first"),
        )
        .reset_index()
    )

    summary["escape_vulnerability_score"] = summary[
        "escape_vulnerability_score"
    ].round(4)

    summary["escape_resilience_score"] = (
        100 * (1 - summary["escape_vulnerability_score"])
    ).round(2)

    summary["escape_resilience_class"] = pd.cut(
        summary["escape_resilience_score"],
        bins=[-1, 50, 70, 100],
        labels=[
            "Low_escape_resilience",
            "Moderate_escape_resilience",
            "High_escape_resilience",
        ],
    )

    summary = summary.sort_values(
        by=["escape_resilience_score", "unseen_gene_priority_score_v1"],
        ascending=[False, False],
    )

    summary.to_csv(OUTPUT_SUMMARY, index=False)

    total_guides = len(summary)
    total_mutants = len(mutants_df)
    high_resilience_guides = (
        summary["escape_resilience_class"].astype(str) == "High_escape_resilience"
    ).sum()

    top_row = summary.iloc[0]

    with open(OUTPUT_TXT, "w") as f:
        f.write("vanA Single-Mutation Escape Simulation Summary\n")
        f.write("============================================\n\n")
        f.write(f"Total vanA guides analyzed: {total_guides}\n")
        f.write(f"Total single-mutant target sequences simulated: {total_mutants}\n")
        f.write(f"High escape-resilience guides: {high_resilience_guides}\n\n")

        f.write("Top escape-resilient guide:\n")
        f.write(f"- Spacer: {top_row['original_spacer']}\n")
        f.write(f"- PAM: {top_row['pam']}\n")
        f.write(f"- Escape resilience score: {top_row['escape_resilience_score']}\n")
        f.write(f"- Escape resilience class: {top_row['escape_resilience_class']}\n")
        f.write(
            f"- First-pass unseen-gene score: "
            f"{top_row['unseen_gene_priority_score_v1']}\n"
        )
        f.write(
            f"- Recommendation v1: "
            f"{top_row['unseen_gene_recommendation_v1']}\n"
        )

    print(f"Wrote: {OUTPUT_MUTANTS}")
    print(f"Wrote: {OUTPUT_SUMMARY}")
    print(f"Wrote: {OUTPUT_TXT}")
    print()
    print("Escape simulation completed.")
    print(f"Total guides analyzed: {total_guides}")
    print(f"Total single-mutant sequences simulated: {total_mutants}")


if __name__ == "__main__":
    main()