#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

INPUT_GUIDE = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/"
    / "requirement2_guide_level_functional_offtarget_burden.csv"
)

INPUT_HITS = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/"
    / "requirement2_expanded_offtarget_hits_functionally_annotated.csv"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2"
)

OUT_FINAL_GUIDES = OUT_DIR / "requirement2_final_guide_prioritization_table.csv"
OUT_TOP_GUIDES = OUT_DIR / "requirement2_top_guides_by_gene.csv"
OUT_DEPRIORITIZED = OUT_DIR / "requirement2_deprioritized_guides.csv"
OUT_SUMMARY = OUT_DIR / "requirement2_completion_summary.txt"

OUT_DIR.mkdir(parents=True, exist_ok=True)

RECOMMENDATION_ORDER = {
    "Strong_candidate_low_functional_offtarget_burden": 1,
    "Acceptable_candidate_moderate_functional_offtarget_burden": 2,
    "Caution_candidate_elevated_functional_offtarget_burden": 3,
    "Avoid_or_deprioritize_high_functional_offtarget_burden": 4,
}


def clean_recommendation(x):
    mapping = {
        "Strong_candidate_low_functional_offtarget_burden": "Strong low-burden candidate",
        "Acceptable_candidate_moderate_functional_offtarget_burden": "Acceptable moderate-burden candidate",
        "Caution_candidate_elevated_functional_offtarget_burden": "Caution elevated-burden candidate",
        "Avoid_or_deprioritize_high_functional_offtarget_burden": "Avoid/deprioritize high-burden candidate",
    }
    return mapping.get(str(x), str(x))


def infer_gene_col(df):
    for col in ["target_gene", "gene", "amr_gene", "target"]:
        if col in df.columns:
            return col
    raise ValueError("Could not find target gene column.")


def infer_guide_col(df):
    for col in ["guide_id", "selected_guide_id", "guide_name", "query_id", "guide"]:
        if col in df.columns:
            return col
    raise ValueError("Could not find guide ID column.")


def main():
    print("Requirement 2 final guide summary")
    print("=================================")

    if not INPUT_GUIDE.exists():
        raise FileNotFoundError(f"Missing guide-level input file: {INPUT_GUIDE}")

    if not INPUT_HITS.exists():
        raise FileNotFoundError(f"Missing annotated hit-level input file: {INPUT_HITS}")

    guide_df = pd.read_csv(INPUT_GUIDE)
    hit_df = pd.read_csv(INPUT_HITS)

    print(f"Guide-level input shape: {guide_df.shape}")
    print(f"Hit-level annotation input shape: {hit_df.shape}")

    gene_col = infer_gene_col(guide_df)
    guide_col = infer_guide_col(guide_df)

    guide_df["recommendation_order"] = (
        guide_df["functional_annotation_recommendation"]
        .map(RECOMMENDATION_ORDER)
        .fillna(99)
    )

    guide_df["recommendation_clean"] = guide_df[
        "functional_annotation_recommendation"
    ].apply(clean_recommendation)

    guide_df = guide_df.sort_values(
        by=[
            gene_col,
            "recommendation_order",
            "total_functional_burden",
            "high_interest_hits",
            "total_offtarget_hits",
        ],
        ascending=True,
    ).copy()

    guide_df["requirement2_rank_within_gene"] = (
        guide_df.groupby(gene_col).cumcount() + 1
    )

    def status_flag(row):
        rec = row["functional_annotation_recommendation"]

        if rec == "Strong_candidate_low_functional_offtarget_burden":
            return "PRIORITIZE"
        elif rec == "Acceptable_candidate_moderate_functional_offtarget_burden":
            return "KEEP_AS_BACKUP"
        elif rec == "Caution_candidate_elevated_functional_offtarget_burden":
            return "USE_WITH_CAUTION"
        else:
            return "DEPRIORITIZE"

    guide_df["requirement2_decision"] = guide_df.apply(status_flag, axis=1)

    top_guides = guide_df[guide_df["requirement2_rank_within_gene"] <= 2].copy()

    deprioritized = guide_df[
        guide_df["functional_annotation_recommendation"]
        == "Avoid_or_deprioritize_high_functional_offtarget_burden"
    ].copy()

    guide_df.to_csv(OUT_FINAL_GUIDES, index=False)
    top_guides.to_csv(OUT_TOP_GUIDES, index=False)
    deprioritized.to_csv(OUT_DEPRIORITIZED, index=False)

    total_hits = len(hit_df)
    accessions = hit_df["accession"].nunique() if "accession" in hit_df.columns else "NA"
    gff_files_used = (
        hit_df["overlap_gff_file"].nunique()
        if "overlap_gff_file" in hit_df.columns
        else "NA"
    )
    missing_gff = (
        (hit_df["overlap_gff_file"].fillna("") == "").sum()
        if "overlap_gff_file" in hit_df.columns
        else "NA"
    )

    genomic_context_counts = hit_df["genomic_context"].value_counts(dropna=False)
    functional_class_counts = hit_df["functional_class_primary"].value_counts(dropna=False)
    recommendation_counts = guide_df[
        "functional_annotation_recommendation"
    ].value_counts(dropna=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 2 Completion Summary\n")
        f.write("================================\n\n")

        f.write("Final status\n")
        f.write("------------\n")
        f.write(
            "Requirement 2 first-pass status: COMPLETE using GFF/GFF3 + "
            "keyword-derived functional annotation.\n"
        )
        f.write(
            "Database-backed CARD/VFDB/DEG/mobile-element confirmation remains pending.\n\n"
        )

        f.write("What was completed\n")
        f.write("------------------\n")
        f.write("- Annotated the expanded Cas-OFFinder hit-level table using matched GFF files.\n")
        f.write("- Assigned each candidate off-target hit a genomic context.\n")
        f.write("- Assigned each hit a first-pass functional class.\n")
        f.write(
            "- Functional classes included AMR/resistance-related, virulence-related, "
            "essential-like, mobile/plasmid-related, coding-other, and annotated non-CDS-other.\n"
        )
        f.write("- Aggregated hit-level annotations into guide-level functional burden scores.\n")
        f.write("- Ranked guides and assigned prioritization decisions for FOTR-CRISPR selection.\n\n")

        f.write("Final metrics\n")
        f.write("-------------\n")
        f.write(f"Annotated off-target hits: {total_hits}\n")
        f.write(f"Genome accessions represented: {accessions}\n")
        f.write(f"GFF files used: {gff_files_used}\n")
        f.write(f"Rows without overlap GFF file: {missing_gff}\n")
        f.write(f"Guide-level rows: {len(guide_df)}\n")
        f.write(f"Top guide rows exported: {len(top_guides)}\n")
        f.write(f"Deprioritized guide rows exported: {len(deprioritized)}\n\n")

        f.write("Genomic context counts\n")
        f.write("----------------------\n")
        f.write(str(genomic_context_counts))
        f.write("\n\n")

        f.write("Primary functional class counts\n")
        f.write("-------------------------------\n")
        f.write(str(functional_class_counts))
        f.write("\n\n")

        f.write("Guide-level recommendation counts\n")
        f.write("---------------------------------\n")
        f.write(str(recommendation_counts))
        f.write("\n\n")

        f.write("Top guides by gene after Requirement 2\n")
        f.write("-------------------------------------\n")

        top_cols = [
            gene_col,
            guide_col,
            "requirement2_rank_within_gene",
            "total_offtarget_hits",
            "amr_hits",
            "virulence_hits",
            "essential_like_hits",
            "mobile_or_plasmid_hits",
            "high_interest_hits",
            "total_functional_burden",
            "recommendation_clean",
            "requirement2_decision",
        ]

        available_top_cols = [c for c in top_cols if c in top_guides.columns]
        f.write(top_guides[available_top_cols].to_string(index=False))
        f.write("\n\n")

        f.write("Manuscript-safe wording\n")
        f.write("-----------------------\n")
        f.write(
            "We extended the expanded Cas-OFFinder off-target table with matched "
            "GFF-based functional annotation across the 36-accession genome panel. "
            "Each of the 1,997 candidate off-target hits was assigned a genomic context "
            "and first-pass functional class based on overlapping genome features and "
            "keyword-derived biological categories. The annotated hit table identified "
            "1,698 coding-other hits, 104 essential-like hits, 77 AMR/resistance-related "
            "hits, 72 virulence-related hits, 30 annotated non-CDS-other hits, and "
            "16 mobile/plasmid-related hits. These hit-level annotations were aggregated "
            "into a guide-level functional off-target burden table to support "
            "FOTR-CRISPR candidate prioritization. This remains a computational "
            "first-pass annotation and will be strengthened using database-backed "
            "CARD, VFDB, DEG, and mobile-element mapping.\n"
        )

        f.write("\nOutput files\n")
        f.write("------------\n")
        f.write(f"Final guide prioritization table: {OUT_FINAL_GUIDES}\n")
        f.write(f"Top guides by gene: {OUT_TOP_GUIDES}\n")
        f.write(f"Deprioritized guides: {OUT_DEPRIORITIZED}\n")
        f.write(f"Completion summary: {OUT_SUMMARY}\n")

    print("\nWrote files:")
    print(f"- {OUT_FINAL_GUIDES}")
    print(f"- {OUT_TOP_GUIDES}")
    print(f"- {OUT_DEPRIORITIZED}")
    print(f"- {OUT_SUMMARY}")

    print("\nTop guides by gene:")
    print(
        top_guides[
            [
                gene_col,
                guide_col,
                "requirement2_rank_within_gene",
                "total_offtarget_hits",
                "high_interest_hits",
                "total_functional_burden",
                "recommendation_clean",
                "requirement2_decision",
            ]
        ].to_string(index=False)
    )

    print("\nRecommendation counts:")
    print(recommendation_counts)

    print("\nDone.")


if __name__ == "__main__":
    main()
