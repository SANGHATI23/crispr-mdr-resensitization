#!/usr/bin/env python3

"""
Requirement 4 Final Integrated Conservation Table

Purpose:
Create a final manuscript-ready all-20 guide table combining:

- Requirement 3 final integrated mobile-context evidence
- Requirement 4 off-target conservation across 36 accessions

This produces the final conserved off-target prioritization table.
"""

from pathlib import Path
import pandas as pd
import numpy as np


PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

REQ3_TABLE = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/mobile_context_mapping/"
    "final_requirement3_complete/"
    "requirement3_final_integrated_mobile_element_mapping_ALL20.csv"
)

REQ4_TABLE = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/"
    "requirement4_guide_level_offtarget_conservation_ALL20.csv"
)

OUT_DIR = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/"
    "final_requirement4_complete"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


def final_conservation_decision(row):
    """
    Conservative final decision after adding Requirement 4.
    Lower score is better.
    """

    total_hits = row.get("total_offtarget_hit_rows", 0)
    req3_score = row.get("requirement3_integrated_evidence_score", 0)
    req4_score = row.get("requirement4_conservation_integrated_score", 0)
    conserved_burden = row.get("conserved_functional_burden_score", 0)
    frac = row.get("fraction_accessions_with_any_hit", 0)

    if total_hits == 0 and req3_score == 0 and req4_score == 0:
        return "FINAL_REQ4_PRIORITIZE_ZERO_HIT_NO_FUNCTIONAL_MOBILE_OR_CONSERVED_RISK"

    if req3_score <= 1 and req4_score <= 1 and conserved_burden <= 25:
        return "FINAL_REQ4_PRIORITIZE_LOW_FUNCTIONAL_MOBILE_AND_CONSERVATION_RISK"

    if req3_score <= 3 and req4_score <= 3 and conserved_burden <= 75:
        return "FINAL_REQ4_BACKUP_OR_USE_WITH_LOW_TO_MODERATE_CONSERVATION_RISK"

    if req3_score <= 5 and req4_score <= 5 and frac < 0.50:
        return "FINAL_REQ4_CAUTION_MODERATE_TO_HIGH_CONSERVED_OFFTARGET_RISK"

    return "FINAL_REQ4_DEPRIORITIZE_HIGH_FUNCTIONAL_MOBILE_OR_CONSERVED_OFFTARGET_RISK"


def decision_order(decision):
    order = {
        "FINAL_REQ4_PRIORITIZE_ZERO_HIT_NO_FUNCTIONAL_MOBILE_OR_CONSERVED_RISK": 1,
        "FINAL_REQ4_PRIORITIZE_LOW_FUNCTIONAL_MOBILE_AND_CONSERVATION_RISK": 2,
        "FINAL_REQ4_BACKUP_OR_USE_WITH_LOW_TO_MODERATE_CONSERVATION_RISK": 3,
        "FINAL_REQ4_CAUTION_MODERATE_TO_HIGH_CONSERVED_OFFTARGET_RISK": 4,
        "FINAL_REQ4_DEPRIORITIZE_HIGH_FUNCTIONAL_MOBILE_OR_CONSERVED_OFFTARGET_RISK": 5,
    }
    return order.get(decision, 99)


def main():
    print("Creating final Requirement 4 integrated conservation table...")

    if not REQ3_TABLE.exists():
        raise FileNotFoundError(f"Missing Requirement 3 table: {REQ3_TABLE}")

    if not REQ4_TABLE.exists():
        raise FileNotFoundError(f"Missing Requirement 4 table: {REQ4_TABLE}")

    req3 = pd.read_csv(REQ3_TABLE)
    req4 = pd.read_csv(REQ4_TABLE)

    print(f"Loaded Requirement 3 table: {req3.shape}")
    print(f"Loaded Requirement 4 table: {req4.shape}")

    key_cols = ["target_gene", "guide_id"]

    keep_req3_cols = [
        "target_gene",
        "guide_id",
        "total_offtarget_hits",
        "db_refined_total_burden",
        "db_refined_high_confidence_functional_hits",
        "requirement3_mobile_context_burden",
        "direct_plasmid_contig_hits",
        "direct_sccmec_relevant_hits",
        "direct_any_mobile_element_marker_hits",
        "requirement3_integrated_evidence_score",
        "requirement3_integrated_confidence_label",
        "requirement3_final_integrated_decision",
        "requirement3_final_integrated_rank_within_gene",
    ]

    keep_req3_cols = [c for c in keep_req3_cols if c in req3.columns]

    keep_req4_cols = [
        "target_gene",
        "guide_id",
        "spacer",
        "total_offtarget_hit_rows",
        "accessions_with_any_hit",
        "total_accessions_screened",
        "fraction_accessions_with_any_hit",
        "guide_conservation_label",
        "unique_offtarget_site_keys",
        "high_confidence_functional_hit_rows",
        "total_db_refined_burden",
        "conserved_functional_burden_score",
        "most_conserved_offtarget_site_key",
        "max_accessions_for_single_offtarget_site",
        "max_single_site_conservation_fraction",
        "most_conserved_site_label",
        "requirement4_conservation_integrated_score",
        "requirement4_final_conservation_decision",
    ]

    keep_req4_cols = [c for c in keep_req4_cols if c in req4.columns]

    final = req4[keep_req4_cols].merge(
        req3[keep_req3_cols],
        on=key_cols,
        how="left",
    )

    # Fill numeric missing values safely.
    numeric_cols = final.select_dtypes(include=[np.number]).columns
    final[numeric_cols] = final[numeric_cols].fillna(0)

    # Final combined score.
    final["final_req2_req3_req4_integrated_score"] = (
        final.get("requirement3_integrated_evidence_score", 0).astype(float)
        + final.get("requirement4_conservation_integrated_score", 0).astype(float)
    )

    final["final_req4_integrated_decision"] = final.apply(
        final_conservation_decision,
        axis=1,
    )

    final["final_req4_decision_order"] = final[
        "final_req4_integrated_decision"
    ].apply(decision_order)

    final = final.sort_values(
        [
            "target_gene",
            "final_req4_decision_order",
            "final_req2_req3_req4_integrated_score",
            "conserved_functional_burden_score",
            "total_offtarget_hit_rows",
        ],
        ascending=[True, True, True, True, True],
    )

    final["final_req4_rank_within_gene"] = (
        final.groupby("target_gene").cumcount() + 1
    )

    final_table = OUT_DIR / "requirement4_final_integrated_conserved_offtarget_table_ALL20.csv"
    final.to_csv(final_table, index=False)

    top_by_gene = final.groupby("target_gene").head(2).copy()
    top_table = OUT_DIR / "requirement4_final_top_guides_by_gene_conservation_integrated.csv"
    top_by_gene.to_csv(top_table, index=False)

    deprioritized = final[
        final["final_req4_integrated_decision"].str.contains("DEPRIORITIZE", na=False)
    ].copy()

    deprioritized_table = OUT_DIR / "requirement4_final_deprioritized_conserved_offtarget_guides.csv"
    deprioritized.to_csv(deprioritized_table, index=False)

    summary_path = OUT_DIR / "requirement4_final_integrated_completion_summary.txt"

    with open(summary_path, "w") as f:
        f.write("Requirement 4 Final Integrated Conservation Summary\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Input Requirement 3 table: {REQ3_TABLE}\n")
        f.write(f"Input Requirement 4 table: {REQ4_TABLE}\n\n")

        f.write(f"Final guide rows retained: {len(final)}\n")
        f.write(f"Unique guides retained: {final['guide_id'].nunique()}\n")
        f.write(f"Target genes represented: {final['target_gene'].nunique()}\n\n")

        f.write("Final Requirement 4 decision counts:\n")
        for decision, count in final["final_req4_integrated_decision"].value_counts().items():
            f.write(f"- {decision}: {count}\n")

        f.write("\nTop guides by gene:\n")
        cols = [
            "target_gene",
            "final_req4_rank_within_gene",
            "guide_id",
            "total_offtarget_hit_rows",
            "accessions_with_any_hit",
            "fraction_accessions_with_any_hit",
            "guide_conservation_label",
            "conserved_functional_burden_score",
            "final_req4_integrated_decision",
        ]
        cols = [c for c in cols if c in top_by_gene.columns]

        for _, row in top_by_gene[cols].iterrows():
            f.write(
                f"- {row['target_gene']} rank {row['final_req4_rank_within_gene']}: "
                f"{row['guide_id']} | hits={row['total_offtarget_hit_rows']} | "
                f"accessions_with_hit={row['accessions_with_any_hit']} | "
                f"conservation_fraction={row['fraction_accessions_with_any_hit']:.3f} | "
                f"decision={row['final_req4_integrated_decision']}\n"
            )

        f.write("\nOutput files:\n")
        f.write(f"- {final_table}\n")
        f.write(f"- {top_table}\n")
        f.write(f"- {deprioritized_table}\n")
        f.write(f"- {summary_path}\n\n")

        f.write("Interpretation:\n")
        f.write(
            "This final Requirement 4 table integrates functional off-target burden, "
            "mobile-genetic-context evidence, and conserved off-target signal across "
            "the expanded 36-accession genome panel. These results remain computational "
            "predictions and should not be described as experimental off-target validation.\n"
        )

    print(f"Wrote final integrated Requirement 4 table: {final_table}")
    print(f"Wrote final top guides by gene: {top_table}")
    print(f"Wrote final deprioritized table: {deprioritized_table}")
    print(f"Wrote final Requirement 4 summary: {summary_path}")

    print("\nFinal top guides by gene:")
    preview_cols = [
        "target_gene",
        "final_req4_rank_within_gene",
        "guide_id",
        "total_offtarget_hit_rows",
        "accessions_with_any_hit",
        "fraction_accessions_with_any_hit",
        "guide_conservation_label",
        "conserved_functional_burden_score",
        "final_req4_integrated_decision",
    ]
    preview_cols = [c for c in preview_cols if c in top_by_gene.columns]
    print(top_by_gene[preview_cols].to_string(index=False))

    print("\nRequirement 4 final integration completed successfully.")


if __name__ == "__main__":
    main()
