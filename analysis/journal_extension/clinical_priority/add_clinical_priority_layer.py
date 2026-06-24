#!/usr/bin/env python3

"""
add_clinical_priority_layer.py

Purpose:
Add a clinical-priority interpretation layer to Pareto-ranked vanA guides.

This module does not claim clinical efficacy. It assigns transparent,
target-side clinical-priority metadata based on AMR gene category, drug class,
mobility context, and treatment relevance.

Inputs:
    results_journal_extension/pareto_ranking/vanA_pareto_ranked_guides.csv
    results_journal_extension/unseen_genes/unseen_amr_gene_panel.csv

Outputs:
    results_journal_extension/clinical_priority/vanA_clinical_priority_guides.csv
    results_journal_extension/clinical_priority/vanA_clinical_priority_top_candidates.csv
    results_journal_extension/clinical_priority/vanA_clinical_priority_summary.txt
"""

from pathlib import Path
import pandas as pd


PARETO_INPUT = Path(
    "results_journal_extension/pareto_ranking/vanA_pareto_ranked_guides.csv"
)

PANEL_INPUT = Path(
    "results_journal_extension/unseen_genes/unseen_amr_gene_panel.csv"
)

OUTPUT_DIR = Path("results_journal_extension/clinical_priority")
OUTPUT_ALL = OUTPUT_DIR / "vanA_clinical_priority_guides.csv"
OUTPUT_TOP = OUTPUT_DIR / "vanA_clinical_priority_top_candidates.csv"
OUTPUT_TXT = OUTPUT_DIR / "vanA_clinical_priority_summary.txt"


def normalize_gene(gene: str) -> str:
    return str(gene).strip().lower()


def clinical_priority_score(gene: str) -> float:
    """
    First-pass target-side clinical priority.

    This is gene-level clinical relevance, not experimental efficacy.
    """
    gene = normalize_gene(gene)

    priority_map = {
        "vana": 1.00,
        "aac6ib": 0.95,
        "qnrs": 0.85,
        "tetm": 0.75,
        "ermb": 0.75,
    }

    return priority_map.get(gene, 0.50)


def treatment_relevance_score(gene: str) -> float:
    """
    First-pass treatment relevance score.
    """
    gene = normalize_gene(gene)

    treatment_map = {
        "vana": 1.00,
        "aac6ib": 0.90,
        "qnrs": 0.85,
        "tetm": 0.70,
        "ermb": 0.70,
    }

    return treatment_map.get(gene, 0.50)


def mobility_context_score(expected_target_context: str) -> float:
    """
    Score whether the target is expected to occur in mobile AMR context.

    This is broad target-context metadata, not confirmed plasmid/genome-wide
    annotation.
    """
    context = str(expected_target_context).lower()

    if "mobile" in context:
        return 1.00
    if "plasmid" in context:
        return 0.95
    if "integron" in context:
        return 0.95
    if "conjugative" in context:
        return 0.90
    if "cluster" in context:
        return 0.85

    return 0.60


def classify_clinical_priority(score: float) -> str:
    if score >= 90:
        return "Very_high_clinical_priority"
    if score >= 75:
        return "High_clinical_priority"
    if score >= 60:
        return "Moderate_clinical_priority"
    return "Lower_clinical_priority"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not PARETO_INPUT.exists():
        raise FileNotFoundError(f"Missing Pareto input: {PARETO_INPUT}")

    if not PANEL_INPUT.exists():
        raise FileNotFoundError(f"Missing unseen-gene panel input: {PANEL_INPUT}")

    pareto = pd.read_csv(PARETO_INPUT)
    panel = pd.read_csv(PANEL_INPUT)

    required_pareto_cols = [
        "gene",
        "spacer",
        "pam",
        "position_1based",
        "target_region_bin",
        "unseen_gene_priority_score_v1",
        "escape_resilience_score",
        "escape_resilience_class",
        "pareto_support_score",
        "is_pareto_front",
        "pareto_recommendation",
    ]

    missing_pareto_cols = [
        col for col in required_pareto_cols if col not in pareto.columns
    ]

    if missing_pareto_cols:
        raise ValueError(f"Missing Pareto columns: {missing_pareto_cols}")

    required_panel_cols = [
        "gene",
        "resistance_category",
        "primary_drug_class",
        "clinical_priority",
        "expected_target_context",
    ]

    missing_panel_cols = [
        col for col in required_panel_cols if col not in panel.columns
    ]

    if missing_panel_cols:
        raise ValueError(f"Missing panel columns: {missing_panel_cols}")

    pareto["gene_normalized"] = pareto["gene"].apply(normalize_gene)
    panel["gene_normalized"] = panel["gene"].apply(normalize_gene)

    merged = pareto.merge(
        panel[
            [
                "gene_normalized",
                "resistance_category",
                "primary_drug_class",
                "clinical_priority",
                "expected_target_context",
            ]
        ],
        on="gene_normalized",
        how="left",
    )

    if merged["resistance_category"].isna().any():
        missing_count = merged["resistance_category"].isna().sum()
        raise ValueError(
            f"{missing_count} rows did not receive clinical panel metadata."
        )

    merged["clinical_priority_component"] = merged["gene"].apply(
        clinical_priority_score
    )

    merged["treatment_relevance_component"] = merged["gene"].apply(
        treatment_relevance_score
    )

    merged["mobility_context_component"] = merged["expected_target_context"].apply(
        mobility_context_score
    )

    merged["clinical_priority_score"] = (
        100
        * (
            0.45 * merged["clinical_priority_component"]
            + 0.35 * merged["treatment_relevance_component"]
            + 0.20 * merged["mobility_context_component"]
        )
    ).round(2)

    merged["clinical_priority_class"] = merged["clinical_priority_score"].apply(
        classify_clinical_priority
    )

    merged["clinical_pareto_support_score"] = (
        0.60 * merged["pareto_support_score"]
        + 0.40 * merged["clinical_priority_score"]
    ).round(2)

    merged["final_clinical_pareto_label"] = "Non_clinical_pareto_candidate"

    merged.loc[
        merged["is_pareto_front"].astype(bool),
        "final_clinical_pareto_label",
    ] = "Clinical_pareto_candidate"

    merged.loc[
        (merged["is_pareto_front"].astype(bool))
        & (merged["pareto_recommendation"] == "High_priority_pareto_candidate")
        & (
            merged["clinical_priority_class"].isin(
                ["Very_high_clinical_priority", "High_clinical_priority"]
            )
        ),
        "final_clinical_pareto_label",
    ] = "High_priority_clinical_pareto_candidate"

    label_priority = {
        "High_priority_clinical_pareto_candidate": 1,
        "Clinical_pareto_candidate": 2,
        "Non_clinical_pareto_candidate": 3,
    }

    merged["label_sort_order"] = merged["final_clinical_pareto_label"].map(
        label_priority
    )

    merged = merged.sort_values(
        by=[
            "label_sort_order",
            "clinical_pareto_support_score",
            "pareto_support_score",
            "clinical_priority_score",
        ],
        ascending=[True, False, False, False],
    ).reset_index(drop=True)

    merged["clinical_pareto_rank"] = range(1, len(merged) + 1)

    output_cols = [
        "clinical_pareto_rank",
        "gene",
        "spacer",
        "pam",
        "position_1based",
        "target_region_bin",
        "resistance_category",
        "primary_drug_class",
        "clinical_priority",
        "expected_target_context",
        "unseen_gene_priority_score_v1",
        "escape_resilience_score",
        "escape_resilience_class",
        "pareto_support_score",
        "pareto_recommendation",
        "clinical_priority_score",
        "clinical_priority_class",
        "clinical_pareto_support_score",
        "final_clinical_pareto_label",
    ]

    merged[output_cols].to_csv(OUTPUT_ALL, index=False)

    top_candidates = merged[
        merged["final_clinical_pareto_label"]
        == "High_priority_clinical_pareto_candidate"
    ].copy()

    top_candidates[output_cols].to_csv(OUTPUT_TOP, index=False)

    total_guides = len(merged)
    pareto_count = merged["is_pareto_front"].astype(bool).sum()
    top_count = len(top_candidates)

    top = merged.iloc[0]

    with open(OUTPUT_TXT, "w") as f:
        f.write("vanA Clinical-Priority Layer Summary\n")
        f.write("====================================\n\n")

        f.write(f"Total vanA guides analyzed: {total_guides}\n")
        f.write(f"Pareto-front candidates carried forward: {pareto_count}\n")
        f.write(f"High-priority clinical Pareto candidates: {top_count}\n\n")

        f.write("Clinical-priority components:\n")
        f.write("- clinical_priority_component\n")
        f.write("- treatment_relevance_component\n")
        f.write("- mobility_context_component\n\n")

        f.write("Top clinical-priority Pareto guide:\n")
        f.write(f"- Spacer: {top['spacer']}\n")
        f.write(f"- PAM: {top['pam']}\n")
        f.write(f"- Position: {top['position_1based']}\n")
        f.write(f"- Target region: {top['target_region_bin']}\n")
        f.write(f"- Resistance category: {top['resistance_category']}\n")
        f.write(f"- Primary drug class: {top['primary_drug_class']}\n")
        f.write(f"- Expected target context: {top['expected_target_context']}\n")
        f.write(f"- Pareto support score: {top['pareto_support_score']}\n")
        f.write(f"- Clinical priority score: {top['clinical_priority_score']}\n")
        f.write(
            f"- Clinical Pareto support score: "
            f"{top['clinical_pareto_support_score']}\n"
        )
        f.write(f"- Final label: {top['final_clinical_pareto_label']}\n")

        f.write("\nBoundary statement:\n")
        f.write(
            "This layer is a target-side clinical-priority interpretation "
            "module. It does not claim experimental cleavage, antimicrobial "
            "resensitization, clinical efficacy, or patient-level utility.\n"
        )

    print(f"Wrote: {OUTPUT_ALL}")
    print(f"Wrote: {OUTPUT_TOP}")
    print(f"Wrote: {OUTPUT_TXT}")
    print()
    print("Clinical-priority layer completed.")
    print(f"Total vanA guides analyzed: {total_guides}")
    print(f"Pareto-front candidates carried forward: {pareto_count}")
    print(f"High-priority clinical Pareto candidates: {top_count}")


if __name__ == "__main__":
    main()