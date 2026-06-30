#!/usr/bin/env python3

"""
fotr_journal_cli.py

Purpose:
Simple command-line interface for inspecting journal-extension FOTR-CRISPR
outputs for prospective unseen-gene validation.

Current supported gene:
    vanA

Usage:
    python analysis/journal_extension/interactive_cli/fotr_journal_cli.py --summary
    python analysis/journal_extension/interactive_cli/fotr_journal_cli.py --top 10
    python analysis/journal_extension/interactive_cli/fotr_journal_cli.py --pareto
    python analysis/journal_extension/interactive_cli/fotr_journal_cli.py --clinical
"""

from pathlib import Path
import argparse
import pandas as pd


GUIDE_INPUT = Path(
    "results_journal_extension/unseen_genes/"
    "unseen_gene_candidate_guides_ranked_v1.csv"
)

ESCAPE_INPUT = Path(
    "results_journal_extension/escape_simulation/"
    "vanA_escape_summary_by_guide.csv"
)

PARETO_INPUT = Path(
    "results_journal_extension/pareto_ranking/"
    "vanA_pareto_ranked_guides.csv"
)

CLINICAL_INPUT = Path(
    "results_journal_extension/clinical_priority/"
    "vanA_clinical_priority_guides.csv"
)


def load_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {label} file: {path}\n"
            f"Run the corresponding journal-extension module first."
        )
    return pd.read_csv(path)


def print_summary() -> None:
    guides = load_csv(GUIDE_INPUT, "ranked unseen-gene guides")
    escape = load_csv(ESCAPE_INPUT, "escape simulation summary")
    pareto = load_csv(PARETO_INPUT, "Pareto ranking")
    clinical = load_csv(CLINICAL_INPUT, "clinical-priority ranking")

    guides_vana = guides[guides["gene"].astype(str).str.lower() == "vana"].copy()

    total_guides = len(guides_vana)
    total_mutants = int(escape["total_single_mutants"].sum())
    pareto_count = int(pareto["is_pareto_front"].astype(bool).sum())
    high_pareto_count = int(
        (pareto["pareto_recommendation"] == "High_priority_pareto_candidate").sum()
    )
    high_clinical_count = int(
        (
            clinical["final_clinical_pareto_label"]
            == "High_priority_clinical_pareto_candidate"
        ).sum()
    )

    top = clinical.iloc[0]

    print()
    print("FOTR-CRISPR Journal-Extension CLI Summary")
    print("========================================")
    print()
    print(f"Gene: vanA")
    print(f"Total ranked candidate guides: {total_guides}")
    print(f"Total simulated single-mutant targets: {total_mutants}")
    print(f"Pareto-front candidates: {pareto_count}")
    print(f"High-priority Pareto candidates: {high_pareto_count}")
    print(f"High-priority clinical Pareto candidates: {high_clinical_count}")
    print()
    print("Top final clinical-priority Pareto guide")
    print("----------------------------------------")
    print(f"Spacer: {top['spacer']}")
    print(f"PAM: {top['pam']}")
    print(f"Position: {top['position_1based']}")
    print(f"Target region: {top['target_region_bin']}")
    print(f"Resistance category: {top['resistance_category']}")
    print(f"Primary drug class: {top['primary_drug_class']}")
    print(f"Escape resilience score: {top['escape_resilience_score']}")
    print(f"Pareto support score: {top['pareto_support_score']}")
    print(f"Clinical priority score: {top['clinical_priority_score']}")
    print(f"Clinical Pareto support score: {top['clinical_pareto_support_score']}")
    print(f"Final label: {top['final_clinical_pareto_label']}")
    print()


def print_top_guides(n: int) -> None:
    clinical = load_csv(CLINICAL_INPUT, "clinical-priority ranking")

    cols = [
        "clinical_pareto_rank",
        "gene",
        "spacer",
        "pam",
        "position_1based",
        "target_region_bin",
        "escape_resilience_class",
        "pareto_support_score",
        "clinical_priority_score",
        "clinical_pareto_support_score",
        "final_clinical_pareto_label",
    ]

    print()
    print(f"Top {n} Clinical-Priority vanA Guides")
    print("====================================")
    print(clinical[cols].head(n).to_string(index=False))
    print()


def print_pareto_guides() -> None:
    pareto = load_csv(PARETO_INPUT, "Pareto ranking")

    pareto_front = pareto[pareto["is_pareto_front"].astype(bool)].copy()

    cols = [
        "pareto_rank",
        "gene",
        "spacer",
        "pam",
        "position_1based",
        "target_region_bin",
        "unseen_gene_priority_score_v1",
        "escape_resilience_score",
        "pareto_support_score",
        "pareto_recommendation",
    ]

    print()
    print("vanA Pareto-Front Candidates")
    print("===========================")
    print(pareto_front[cols].to_string(index=False))
    print()


def print_clinical_guides() -> None:
    clinical = load_csv(CLINICAL_INPUT, "clinical-priority ranking")

    high = clinical[
        clinical["final_clinical_pareto_label"]
        == "High_priority_clinical_pareto_candidate"
    ].copy()

    cols = [
        "clinical_pareto_rank",
        "gene",
        "spacer",
        "pam",
        "position_1based",
        "target_region_bin",
        "resistance_category",
        "primary_drug_class",
        "clinical_pareto_support_score",
        "final_clinical_pareto_label",
    ]

    print()
    print("High-Priority Clinical Pareto Candidates")
    print("========================================")
    print(high[cols].to_string(index=False))
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect FOTR-CRISPR journal-extension outputs."
    )

    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print overall journal-extension summary.",
    )

    parser.add_argument(
        "--top",
        type=int,
        default=None,
        help="Print top N clinical-priority guides.",
    )

    parser.add_argument(
        "--pareto",
        action="store_true",
        help="Print Pareto-front candidates.",
    )

    parser.add_argument(
        "--clinical",
        action="store_true",
        help="Print high-priority clinical Pareto candidates.",
    )

    args = parser.parse_args()

    if args.summary:
        print_summary()

    if args.top is not None:
        print_top_guides(args.top)

    if args.pareto:
        print_pareto_guides()

    if args.clinical:
        print_clinical_guides()

    if not any([args.summary, args.top is not None, args.pareto, args.clinical]):
        parser.print_help()


if __name__ == "__main__":
    main()