#!/usr/bin/env python3

"""
simulate_single_mutation_escape.py

Purpose:
Simulate first-pass, sequence-aware single-nucleotide escape risk for vanA
candidate CRISPR guides from the prospective unseen-gene validation module.

This script mutates every spacer position to the three alternative nucleotides,
assigns a position-aware and sequence-aware escape-risk weight, and summarizes
guide-level escape resilience.

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
    First-pass positional escape-risk classification.

    Positions 1-7 are treated as highest-risk escape positions.
    Positions 8-12 are treated as moderate-risk positions.
    Positions 13-20 are treated as lower-risk positions.

    This is a computational prioritization heuristic, not experimental proof.
    """
    if 1 <= position_1based <= 7:
        return "High_escape_risk"
    if 8 <= position_1based <= 12:
        return "Moderate_escape_risk"
    return "Lower_escape_risk"


def get_risk_weight(
    risk_class: str,
    original_base: str,
    mutated_base: str,
    spacer: str,
    position_1based: int,
) -> float:
    """
    Sequence-aware escape-risk weight.

    This keeps the positional escape-risk logic but adds guide-specific
    sequence effects so not every guide receives the same score.

    Factors included:
        1. Position-based risk class.
        2. GC-to-AT or AT-to-GC mutation effect.
        3. Local repetitive sequence context.
        4. Homopolymer-like spacer pattern.

    The final weight is capped at 1.0.
    """

    if risk_class == "High_escape_risk":
        weight = 1.0
    elif risk_class == "Moderate_escape_risk":
        weight = 0.5
    elif risk_class == "Lower_escape_risk":
        weight = 0.1
    else:
        weight = 0.0

    original_base = str(original_base).upper()
    mutated_base = str(mutated_base).upper()
    spacer = str(spacer).upper()

    original_is_gc = original_base in {"G", "C"}
    mutated_is_gc = mutated_base in {"G", "C"}

    # GC-to-AT or AT-to-GC changes may alter local binding stability.
    if original_is_gc != mutated_is_gc:
        weight += 0.10

    # Local repeated sequence context is treated as slightly more vulnerable.
    start = max(0, position_1based - 3)
    end = min(len(spacer), position_1based + 2)
    local_window = spacer[start:end]

    if local_window:
        max_base_repeat = max(local_window.count(base) for base in BASES)
        if max_base_repeat >= 4:
            weight += 0.10

    # Homopolymer-like guides are slightly penalized.
    homopolymer_patterns = ["AAAA", "CCCC", "GGGG", "TTTT"]
    if any(pattern in spacer for pattern in homopolymer_patterns):
        weight += 0.05

    return min(weight, 1.0)


def generate_single_mutants(row: pd.Series) -> list[dict]:
    """
    Generate all possible single-nucleotide mutant spacers for one guide.

    A 20-nt spacer produces:
        20 positions x 3 alternative bases = 60 mutant spacers.
    """

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
                    "escape_risk_weight": get_risk_weight(
                        escape_class,
                        original_base,
                        mutated_base,
                        spacer,
                        mutation_position,
                    ),
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

    # Keep vanA only for this first prospective unseen-gene implementation.
    if "gene" in guides.columns:
        guides = guides[guides["gene"].astype(str).str.lower() == "vana"].copy()

    if guides.empty:
        raise ValueError(
            "No vanA guide rows found in the input file. "
            "Check the gene column and input path."
        )

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

    # Relative classification is better here than fixed absolute cutoffs,
    # because this simulation intentionally includes many seed-region mutations
    # for every guide. The top third are therefore labeled high relative
    # escape-resilience within the vanA candidate pool.
    summary["escape_resilience_percentile"] = (
        summary["escape_resilience_score"].rank(pct=True) * 100
    ).round(2)

    summary["escape_resilience_class"] = pd.cut(
        summary["escape_resilience_percentile"],
        bins=[-1, 33.33, 66.67, 100],
        labels=[
            "Lower_escape_resilience",
            "Moderate_escape_resilience",
            "High_escape_resilience",
        ],
    )

    summary = summary.sort_values(
        by=[
            "escape_resilience_score",
            "unseen_gene_priority_score_v1",
            "simple_on_target_proxy",
        ],
        ascending=[False, False, False],
    )

    summary.to_csv(OUTPUT_SUMMARY, index=False)

    total_guides = len(summary)
    total_mutants = len(mutants_df)

    high_resilience_guides = (
        summary["escape_resilience_class"].astype(str) == "High_escape_resilience"
    ).sum()

    moderate_resilience_guides = (
        summary["escape_resilience_class"].astype(str) == "Moderate_escape_resilience"
    ).sum()

    lower_resilience_guides = (
        summary["escape_resilience_class"].astype(str) == "Lower_escape_resilience"
    ).sum()

    top_row = summary.iloc[0]

    with open(OUTPUT_TXT, "w") as f:
        f.write("vanA Single-Mutation Escape Simulation Summary\n")
        f.write("============================================\n\n")
        f.write(f"Total vanA guides analyzed: {total_guides}\n")
        f.write(f"Total single-mutant target sequences simulated: {total_mutants}\n\n")

        f.write("Relative escape-resilience distribution:\n")
        f.write(f"- High_escape_resilience: {high_resilience_guides}\n")
        f.write(f"- Moderate_escape_resilience: {moderate_resilience_guides}\n")
        f.write(f"- Lower_escape_resilience: {lower_resilience_guides}\n\n")

        f.write("Top escape-resilient guide:\n")
        f.write(f"- Spacer: {top_row['original_spacer']}\n")
        f.write(f"- PAM: {top_row['pam']}\n")
        f.write(f"- Escape resilience score: {top_row['escape_resilience_score']}\n")
        f.write(
            f"- Escape resilience percentile: "
            f"{top_row['escape_resilience_percentile']}\n"
        )
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
    print(f"High escape-resilience guides: {high_resilience_guides}")
    print(f"Moderate escape-resilience guides: {moderate_resilience_guides}")
    print(f"Lower escape-resilience guides: {lower_resilience_guides}")


if __name__ == "__main__":
    main()