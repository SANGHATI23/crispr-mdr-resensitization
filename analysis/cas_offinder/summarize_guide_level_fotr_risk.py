import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

ANNOTATED_FILE = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "annotated"
    / "cas_offinder_offtarget_hits_gff_annotated.csv"
)

OUTPUT_DIR = PROJECT_ROOT / "results_cas_offinder" / "guide_level_fotr"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "guide_level_functional_offtarget_risk_summary.csv"
TOP_FILE = OUTPUT_DIR / "final_lowest_risk_guides_for_mic_planning.csv"
SUMMARY_FILE = OUTPUT_DIR / "guide_level_fotr_summary.txt"


def calculate_hit_weight(row):
    """
    Higher weight = more concerning off-target hit.

    Logic:
    - Lower mismatch count is more serious.
    - Coding hits are more important than intergenic hits.
    - AMR/resistance, mobile element, and essential-like hits receive extra penalty.
    """

    mismatch = row.get("mismatch_count")

    if pd.isna(mismatch):
        mismatch_weight = 0.5
    elif mismatch <= 1:
        mismatch_weight = 10.0
    elif mismatch == 2:
        mismatch_weight = 6.0
    elif mismatch == 3:
        mismatch_weight = 3.0
    else:
        mismatch_weight = 1.0

    context = str(row.get("genomic_context", ""))

    if context == "Coding_CDS":
        context_weight = 2.0
    elif context == "Gene_region_nonCDS_or_gene_feature":
        context_weight = 1.5
    elif context == "Noncoding_RNA":
        context_weight = 1.2
    elif context == "Intergenic_or_unannotated":
        context_weight = 0.4
    else:
        context_weight = 1.0

    functional_bonus = 0.0

    if bool(row.get("keyword_amr_like", False)):
        functional_bonus += 5.0

    if bool(row.get("keyword_mobile_element_like", False)):
        functional_bonus += 4.0

    if bool(row.get("keyword_essential_like", False)):
        functional_bonus += 6.0

    return mismatch_weight * context_weight + functional_bonus


def recommendation_label(score, total_hits, high_interest_hits):
    """
    Conservative first-pass label.
    You can tune thresholds later.
    """
    if score <= 8 and high_interest_hits == 0:
        return "Strong_candidate_low_functional_offtarget_burden"
    elif score <= 25 and high_interest_hits <= 1:
        return "Acceptable_candidate_moderate_functional_offtarget_burden"
    elif score <= 60:
        return "Caution_candidate_elevated_functional_offtarget_burden"
    else:
        return "Avoid_or_deprioritize_high_functional_offtarget_burden"


def main():
    print("Creating guide-level FOTR risk summary...")

    if not ANNOTATED_FILE.exists():
        raise FileNotFoundError(f"Missing annotated file: {ANNOTATED_FILE}")

    df = pd.read_csv(ANNOTATED_FILE)
    print(f"Loaded annotated off-target table: {df.shape}")

    required = [
        "final_guide_id",
        "mlcb_gene",
        "final_gene_rank",
        "query_spacer",
        "genome_label",
        "mismatch_count",
        "genomic_context",
        "keyword_amr_like",
        "keyword_mobile_element_like",
        "keyword_essential_like",
        "functional_risk_penalized_score",
        "predicted_high_functional_offtarget_risk_probability",
        "first_pass_functional_offtarget_class",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["hit_risk_weight"] = df.apply(calculate_hit_weight, axis=1)

    df["is_high_interest_hit"] = (
        (df["mismatch_count"] <= 3)
        | (df["keyword_amr_like"] == True)
        | (df["keyword_mobile_element_like"] == True)
        | (df["keyword_essential_like"] == True)
    )

    grouped = []

    for guide_id, sub in df.groupby("final_guide_id"):
        row0 = sub.iloc[0]

        total_hits = len(sub)
        high_interest_hits = int(sub["is_high_interest_hit"].sum())

        summary = {
            "final_guide_id": guide_id,
            "mlcb_gene": row0["mlcb_gene"],
            "original_gene_rank": row0["final_gene_rank"],
            "query_spacer": row0["query_spacer"],
            "original_functional_risk_penalized_score": row0["functional_risk_penalized_score"],
            "original_predicted_high_functional_offtarget_risk_probability": row0[
                "predicted_high_functional_offtarget_risk_probability"
            ],
            "total_offtarget_hits": total_hits,
            "high_interest_hits": high_interest_hits,
            "mismatch_1_hits": int((sub["mismatch_count"] == 1).sum()),
            "mismatch_2_hits": int((sub["mismatch_count"] == 2).sum()),
            "mismatch_3_hits": int((sub["mismatch_count"] == 3).sum()),
            "mismatch_4_hits": int((sub["mismatch_count"] == 4).sum()),
            "coding_cds_hits": int((sub["genomic_context"] == "Coding_CDS").sum()),
            "intergenic_or_unannotated_hits": int(
                (sub["genomic_context"] == "Intergenic_or_unannotated").sum()
            ),
            "amr_or_resistance_keyword_hits": int((sub["keyword_amr_like"] == True).sum()),
            "mobile_element_keyword_hits": int((sub["keyword_mobile_element_like"] == True).sum()),
            "essential_like_keyword_hits": int((sub["keyword_essential_like"] == True).sum()),
            "unique_genomes_with_hits": sub["genome_label"].nunique(),
            "functional_offtarget_burden_score": float(sub["hit_risk_weight"].sum()),
        }

        summary["guide_level_fotr_recommendation"] = recommendation_label(
            summary["functional_offtarget_burden_score"],
            total_hits,
            high_interest_hits,
        )

        grouped.append(summary)

    guide_summary = pd.DataFrame(grouped)

    # Lower burden score is better.
    # Higher original penalized score is better.
    # We create a final rank within each gene.
    guide_summary = guide_summary.sort_values(
        by=[
            "mlcb_gene",
            "functional_offtarget_burden_score",
            "high_interest_hits",
            "total_offtarget_hits",
            "original_functional_risk_penalized_score",
        ],
        ascending=[True, True, True, True, False],
    ).copy()

    guide_summary["final_fotr_rank_within_gene"] = (
        guide_summary.groupby("mlcb_gene").cumcount() + 1
    )

    # Pick top 2 per gene for future MIC planning.
    top_for_mic = guide_summary[guide_summary["final_fotr_rank_within_gene"] <= 2].copy()

    guide_summary.to_csv(OUTPUT_FILE, index=False)
    top_for_mic.to_csv(TOP_FILE, index=False)

    with open(SUMMARY_FILE, "w") as f:
        f.write("Guide-Level FOTR Risk Summary\n")
        f.write("=============================\n\n")
        f.write(f"Input annotated off-target table: {ANNOTATED_FILE}\n")
        f.write(f"Total annotated hits: {df.shape[0]}\n")
        f.write(f"Guides with at least one off-target hit: {guide_summary.shape[0]}\n\n")

        f.write("Guide-level FOTR recommendations:\n")
        rec_counts = (
            guide_summary["guide_level_fotr_recommendation"]
            .value_counts()
            .reset_index()
        )
        rec_counts.columns = ["recommendation", "guide_count"]
        f.write(rec_counts.to_string(index=False))
        f.write("\n\n")

        f.write("Top 2 guides per gene for MIC planning:\n")
        f.write(
            top_for_mic[
                [
                    "mlcb_gene",
                    "final_fotr_rank_within_gene",
                    "final_guide_id",
                    "query_spacer",
                    "functional_offtarget_burden_score",
                    "total_offtarget_hits",
                    "high_interest_hits",
                    "guide_level_fotr_recommendation",
                ]
            ].to_string(index=False)
        )
        f.write("\n\n")

        f.write("Full guide-level table:\n")
        f.write(
            guide_summary[
                [
                    "mlcb_gene",
                    "final_fotr_rank_within_gene",
                    "final_guide_id",
                    "query_spacer",
                    "functional_offtarget_burden_score",
                    "total_offtarget_hits",
                    "high_interest_hits",
                    "mismatch_1_hits",
                    "mismatch_2_hits",
                    "mismatch_3_hits",
                    "mismatch_4_hits",
                    "amr_or_resistance_keyword_hits",
                    "mobile_element_keyword_hits",
                    "essential_like_keyword_hits",
                    "guide_level_fotr_recommendation",
                ]
            ].to_string(index=False)
        )

    print("\nWrote:")
    print(f"- {OUTPUT_FILE}")
    print(f"- {TOP_FILE}")
    print(f"- {SUMMARY_FILE}")

    print("\nTop 2 guides per gene for MIC planning:")
    print(
        top_for_mic[
            [
                "mlcb_gene",
                "final_fotr_rank_within_gene",
                "final_guide_id",
                "query_spacer",
                "functional_offtarget_burden_score",
                "total_offtarget_hits",
                "high_interest_hits",
                "guide_level_fotr_recommendation",
            ]
        ].to_string(index=False)
    )

    print("\nDone.")


if __name__ == "__main__":
    main()