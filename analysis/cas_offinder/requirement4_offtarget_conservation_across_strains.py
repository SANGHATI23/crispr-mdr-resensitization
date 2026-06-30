#!/usr/bin/env python3

"""
Requirement 4: Off-target conservation across many strain genomes

Purpose:
This script calculates whether predicted off-target hits are conserved/repeated
across many genome accessions/strains.

This is the missing Requirement 4 layer:
- not just "does a guide have off-targets?"
- but "are those off-targets repeatedly observed across many strains?"

Inputs:
1. Requirement 2 database-refined hit annotation table
2. Requirement 3 final integrated guide table

Outputs:
1. Hit-level conservation-ready table
2. Off-target-site-level conservation table
3. Guide-level strain conservation summary
4. Requirement 4 completion summary text
"""

from pathlib import Path
import pandas as pd
import numpy as np


PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

REQ2_HIT_TABLE = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/functional_annotation/"
    "final_requirement2_database_refined/"
    "requirement2_database_refined_hit_annotations.csv"
)

REQ3_GUIDE_TABLE = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/mobile_context_mapping/"
    "final_requirement3_complete/"
    "requirement3_final_integrated_mobile_element_mapping_ALL20.csv"
)

OUT_DIR = PROJECT_ROOT / (
    "results_cas_offinder/expanded_panel/offtarget_conservation_requirement4"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


def classify_conservation_fraction(frac):
    """
    Classify conservation level based on fraction of screened accessions.
    """
    if pd.isna(frac) or frac == 0:
        return "No_detected_offtarget_conservation"
    elif frac < 0.10:
        return "Low_conservation_rare_across_panel"
    elif frac < 0.25:
        return "Moderate_conservation_across_panel"
    elif frac < 0.50:
        return "High_conservation_across_panel"
    else:
        return "Very_high_conservation_across_panel"


def conservation_risk_points(frac):
    """
    Convert conservation fraction into simple risk points.
    Higher means more concerning.
    """
    if pd.isna(frac) or frac == 0:
        return 0
    elif frac < 0.10:
        return 1
    elif frac < 0.25:
        return 2
    elif frac < 0.50:
        return 3
    else:
        return 4


def main():
    print("Starting Requirement 4: off-target conservation across strains...")

    if not REQ2_HIT_TABLE.exists():
        raise FileNotFoundError(f"Missing Requirement 2 hit table: {REQ2_HIT_TABLE}")

    if not REQ3_GUIDE_TABLE.exists():
        raise FileNotFoundError(f"Missing Requirement 3 guide table: {REQ3_GUIDE_TABLE}")

    hits = pd.read_csv(REQ2_HIT_TABLE)
    guides = pd.read_csv(REQ3_GUIDE_TABLE)

    print(f"Loaded Requirement 2 hit table: {hits.shape}")
    print(f"Loaded Requirement 3 guide table: {guides.shape}")

    # Confirm required columns using your actual column names.
    required_hit_cols = [
        "target_gene",
        "guide_id",
        "spacer",
        "accession",
        "contig_description",
        "position_1based",
        "strand",
        "off_target_sequence",
        "mismatch_count",
        "genomic_context",
        "db_refined_high_confidence_functional_hit",
        "db_refined_burden_score",
    ]

    missing_cols = [c for c in required_hit_cols if c not in hits.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in hit table: {missing_cols}")

    total_accessions = hits["accession"].nunique()
    total_guides_in_hits = hits["guide_id"].nunique()

    print(f"Unique accessions represented in hit table: {total_accessions}")
    print(f"Guides with at least one hit: {total_guides_in_hits}")

    # ---------------------------------------------------------------------
    # 1. Create hit-level conservation-ready table
    # ---------------------------------------------------------------------
    conservation_ready = hits.copy()

    conservation_ready["offtarget_site_key"] = (
        conservation_ready["guide_id"].astype(str)
        + "|"
        + conservation_ready["off_target_sequence"].astype(str)
        + "|mm"
        + conservation_ready["mismatch_count"].astype(str)
    )

    conservation_ready["guide_accession_key"] = (
        conservation_ready["guide_id"].astype(str)
        + "|"
        + conservation_ready["accession"].astype(str)
    )

    hit_level_out = OUT_DIR / "requirement4_hit_level_conservation_ready_table.csv"
    conservation_ready.to_csv(hit_level_out, index=False)

    print(f"Wrote hit-level conservation-ready table: {hit_level_out}")

    # ---------------------------------------------------------------------
    # 2. Off-target-site-level conservation
    # ---------------------------------------------------------------------
    site_level = (
        conservation_ready
        .groupby(
            [
                "target_gene",
                "guide_id",
                "spacer",
                "offtarget_site_key",
                "off_target_sequence",
                "mismatch_count",
            ],
            dropna=False,
        )
        .agg(
            accessions_with_this_site=("accession", "nunique"),
            total_rows_for_this_site=("accession", "size"),
            min_position_1based=("position_1based", "min"),
            max_position_1based=("position_1based", "max"),
            high_confidence_rows=(
                "db_refined_high_confidence_functional_hit",
                lambda x: int(pd.Series(x).fillna(False).astype(bool).sum()),
            ),
            mean_db_refined_burden=("db_refined_burden_score", "mean"),
            max_db_refined_burden=("db_refined_burden_score", "max"),
        )
        .reset_index()
    )

    site_level["fraction_of_accessions_with_site"] = (
        site_level["accessions_with_this_site"] / total_accessions
    )

    site_level["site_conservation_label"] = site_level[
        "fraction_of_accessions_with_site"
    ].apply(classify_conservation_fraction)

    site_level["site_conservation_risk_points"] = site_level[
        "fraction_of_accessions_with_site"
    ].apply(conservation_risk_points)

    site_level = site_level.sort_values(
        [
            "site_conservation_risk_points",
            "accessions_with_this_site",
            "max_db_refined_burden",
            "high_confidence_rows",
        ],
        ascending=[False, False, False, False],
    )

    site_level_out = OUT_DIR / "requirement4_offtarget_site_level_conservation.csv"
    site_level.to_csv(site_level_out, index=False)

    print(f"Wrote site-level conservation table: {site_level_out}")

    # ---------------------------------------------------------------------
    # 3. Guide-level conservation summary
    # ---------------------------------------------------------------------
    guide_level = (
        conservation_ready
        .groupby(["target_gene", "guide_id", "spacer"], dropna=False)
        .agg(
            total_offtarget_hit_rows=("accession", "size"),
            accessions_with_any_hit=("accession", "nunique"),
            unique_offtarget_site_keys=("offtarget_site_key", "nunique"),
            high_confidence_functional_hit_rows=(
                "db_refined_high_confidence_functional_hit",
                lambda x: int(pd.Series(x).fillna(False).astype(bool).sum()),
            ),
            total_db_refined_burden=("db_refined_burden_score", "sum"),
            mean_db_refined_burden=("db_refined_burden_score", "mean"),
            max_db_refined_burden=("db_refined_burden_score", "max"),
            min_mismatch_count=("mismatch_count", "min"),
            max_mismatch_count=("mismatch_count", "max"),
        )
        .reset_index()
    )

    guide_level["total_accessions_screened"] = total_accessions

    guide_level["fraction_accessions_with_any_hit"] = (
        guide_level["accessions_with_any_hit"] / total_accessions
    )

    guide_level["guide_conservation_label"] = guide_level[
        "fraction_accessions_with_any_hit"
    ].apply(classify_conservation_fraction)

    guide_level["guide_conservation_risk_points"] = guide_level[
        "fraction_accessions_with_any_hit"
    ].apply(conservation_risk_points)

    # Conserved functional burden:
    # combines how often the guide has off-targets across strains with
    # its Requirement 2 database-refined burden.
    guide_level["conserved_functional_burden_score"] = (
        guide_level["total_db_refined_burden"]
        * (1 + guide_level["fraction_accessions_with_any_hit"])
    )

    # Add site-level strongest conserved site per guide.
    strongest_site = (
        site_level
        .sort_values(
            [
                "guide_id",
                "site_conservation_risk_points",
                "accessions_with_this_site",
                "max_db_refined_burden",
            ],
            ascending=[True, False, False, False],
        )
        .groupby("guide_id", as_index=False)
        .first()
        [
            [
                "guide_id",
                "offtarget_site_key",
                "accessions_with_this_site",
                "fraction_of_accessions_with_site",
                "site_conservation_label",
                "site_conservation_risk_points",
            ]
        ]
        .rename(
            columns={
                "offtarget_site_key": "most_conserved_offtarget_site_key",
                "accessions_with_this_site": "max_accessions_for_single_offtarget_site",
                "fraction_of_accessions_with_site": "max_single_site_conservation_fraction",
                "site_conservation_label": "most_conserved_site_label",
                "site_conservation_risk_points": "most_conserved_site_risk_points",
            }
        )
    )

    guide_level = guide_level.merge(strongest_site, on="guide_id", how="left")

    # ---------------------------------------------------------------------
    # 4. Restore all 20 guides using Requirement 3 guide table
    # ---------------------------------------------------------------------
    all_guides = guides[["target_gene", "guide_id"]].drop_duplicates()

    guide_level_all20 = all_guides.merge(
        guide_level,
        on=["target_gene", "guide_id"],
        how="left",
    )

    # Fill zero-hit guide values.
    numeric_fill_zero = [
        "total_offtarget_hit_rows",
        "accessions_with_any_hit",
        "unique_offtarget_site_keys",
        "high_confidence_functional_hit_rows",
        "total_db_refined_burden",
        "mean_db_refined_burden",
        "max_db_refined_burden",
        "guide_conservation_risk_points",
        "conserved_functional_burden_score",
        "max_accessions_for_single_offtarget_site",
        "max_single_site_conservation_fraction",
        "most_conserved_site_risk_points",
    ]

    for col in numeric_fill_zero:
        if col in guide_level_all20.columns:
            guide_level_all20[col] = guide_level_all20[col].fillna(0)

    if "total_accessions_screened" in guide_level_all20.columns:
        guide_level_all20["total_accessions_screened"] = guide_level_all20[
            "total_accessions_screened"
        ].fillna(total_accessions)
    else:
        guide_level_all20["total_accessions_screened"] = total_accessions

    guide_level_all20["fraction_accessions_with_any_hit"] = guide_level_all20[
        "fraction_accessions_with_any_hit"
    ].fillna(0)

    guide_level_all20["guide_conservation_label"] = guide_level_all20[
        "guide_conservation_label"
    ].fillna("No_detected_offtarget_conservation")

    guide_level_all20["most_conserved_site_label"] = guide_level_all20[
        "most_conserved_site_label"
    ].fillna("No_detected_offtarget_site")

    guide_level_all20["most_conserved_offtarget_site_key"] = guide_level_all20[
        "most_conserved_offtarget_site_key"
    ].fillna("No_detected_offtarget_site")

    # Add Requirement 3 final integrated fields.
    req3_keep_cols = [
        "target_gene",
        "guide_id",
        "requirement3_integrated_evidence_score",
        "requirement3_integrated_confidence_label",
        "requirement3_final_integrated_decision",
        "requirement3_final_decision_order",
    ]

    req3_keep_cols = [c for c in req3_keep_cols if c in guides.columns]

    guide_level_all20 = guide_level_all20.merge(
        guides[req3_keep_cols].drop_duplicates(),
        on=["target_gene", "guide_id"],
        how="left",
    )

    # Final Requirement 4 decision score.
    guide_level_all20["requirement4_conservation_integrated_score"] = (
        guide_level_all20["guide_conservation_risk_points"].astype(float)
        + guide_level_all20["most_conserved_site_risk_points"].astype(float)
        + np.where(
            guide_level_all20["high_confidence_functional_hit_rows"].astype(float) > 0,
            1,
            0,
        )
    )

    def final_req4_decision(row):
        score = row["requirement4_conservation_integrated_score"]
        frac = row["fraction_accessions_with_any_hit"]

        if row["total_offtarget_hit_rows"] == 0:
            return "REQ4_PRIORITIZE_ZERO_HIT_NO_CONSERVED_OFFTARGET"
        if score <= 1 and frac < 0.10:
            return "REQ4_LOW_CONSERVED_OFFTARGET_RISK"
        if score <= 3:
            return "REQ4_MODERATE_CONSERVED_OFFTARGET_RISK"
        if score <= 5:
            return "REQ4_CAUTION_HIGH_CONSERVED_OFFTARGET_RISK"
        return "REQ4_DEPRIORITIZE_VERY_HIGH_CONSERVED_OFFTARGET_RISK"

    guide_level_all20["requirement4_final_conservation_decision"] = (
        guide_level_all20.apply(final_req4_decision, axis=1)
    )

    guide_level_all20 = guide_level_all20.sort_values(
        [
            "target_gene",
            "requirement4_conservation_integrated_score",
            "conserved_functional_burden_score",
            "total_offtarget_hit_rows",
        ],
        ascending=[True, True, True, True],
    )

    guide_level_out = OUT_DIR / "requirement4_guide_level_offtarget_conservation_ALL20.csv"
    guide_level_all20.to_csv(guide_level_out, index=False)

    print(f"Wrote guide-level conservation summary ALL20: {guide_level_out}")

    # ---------------------------------------------------------------------
    # 5. Top guides by gene after Requirement 4
    # ---------------------------------------------------------------------
    top_by_gene = (
        guide_level_all20
        .sort_values(
            [
                "target_gene",
                "requirement4_conservation_integrated_score",
                "conserved_functional_burden_score",
                "total_offtarget_hit_rows",
            ],
            ascending=[True, True, True, True],
        )
        .groupby("target_gene")
        .head(2)
        .copy()
    )

    top_by_gene["requirement4_rank_within_gene"] = (
        top_by_gene
        .groupby("target_gene")
        .cumcount()
        + 1
    )

    top_by_gene_out = OUT_DIR / "requirement4_top_guides_by_gene_offtarget_conservation.csv"
    top_by_gene.to_csv(top_by_gene_out, index=False)

    print(f"Wrote Requirement 4 top guides by gene: {top_by_gene_out}")

    # ---------------------------------------------------------------------
    # 6. Summary text
    # ---------------------------------------------------------------------
    summary_path = OUT_DIR / "requirement4_offtarget_conservation_completion_summary.txt"

    decision_counts = guide_level_all20[
        "requirement4_final_conservation_decision"
    ].value_counts().to_dict()

    label_counts = guide_level_all20[
        "guide_conservation_label"
    ].value_counts().to_dict()

    with open(summary_path, "w") as f:
        f.write("Requirement 4: Off-target conservation across many strain genomes\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Input Requirement 2 hit table: {REQ2_HIT_TABLE}\n")
        f.write(f"Input Requirement 3 guide table: {REQ3_GUIDE_TABLE}\n\n")
        f.write(f"Total hit rows analyzed: {len(hits)}\n")
        f.write(f"Total accessions screened: {total_accessions}\n")
        f.write(f"Guides retained in final ALL20 table: {guide_level_all20['guide_id'].nunique()}\n")
        f.write(f"Unique off-target site keys: {site_level['offtarget_site_key'].nunique()}\n\n")

        f.write("Guide conservation label counts:\n")
        for k, v in label_counts.items():
            f.write(f"- {k}: {v}\n")

        f.write("\nRequirement 4 final decision counts:\n")
        for k, v in decision_counts.items():
            f.write(f"- {k}: {v}\n")

        f.write("\nOutput files:\n")
        f.write(f"- {hit_level_out}\n")
        f.write(f"- {site_level_out}\n")
        f.write(f"- {guide_level_out}\n")
        f.write(f"- {top_by_gene_out}\n")

        f.write("\nInterpretation:\n")
        f.write(
            "This analysis estimates whether predicted off-target hits are repeatedly "
            "observed across the 36-accession expanded genome panel. It provides a "
            "formal guide-level conserved off-target burden layer. These results remain "
            "computational predictions and do not represent experimental off-target validation.\n"
        )

    print(f"Wrote Requirement 4 completion summary: {summary_path}")

    print("\nRequirement 4 Step 2 completed successfully.")
    print("\nTop guides by gene preview:")
    preview_cols = [
        "target_gene",
        "requirement4_rank_within_gene",
        "guide_id",
        "total_offtarget_hit_rows",
        "accessions_with_any_hit",
        "fraction_accessions_with_any_hit",
        "guide_conservation_label",
        "conserved_functional_burden_score",
        "requirement4_final_conservation_decision",
    ]
    preview_cols = [c for c in preview_cols if c in top_by_gene.columns]
    print(top_by_gene[preview_cols].to_string(index=False))


if __name__ == "__main__":
    main()
