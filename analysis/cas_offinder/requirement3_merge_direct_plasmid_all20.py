#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

COMBINED_ALL20 = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_combined/requirement3_combined_requirement2_requirement3_final_guide_table_ALL20.csv"

DIRECT_PLASMID_GUIDE = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_plasmid_mapping/requirement3_direct_plasmid_guide_level_summary.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_direct_evidence"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL20 = OUT_DIR / "requirement3_all20_with_direct_plasmid_evidence.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_direct_plasmid_all20_summary.txt"


def main():
    print("Merging direct plasmid evidence into all-20 Requirement 3 table...")

    if not COMBINED_ALL20.exists():
        raise FileNotFoundError(f"Missing combined all-20 table: {COMBINED_ALL20}")

    if not DIRECT_PLASMID_GUIDE.exists():
        raise FileNotFoundError(f"Missing direct plasmid guide table: {DIRECT_PLASMID_GUIDE}")

    combined = pd.read_csv(COMBINED_ALL20)
    plasmid = pd.read_csv(DIRECT_PLASMID_GUIDE)

    print(f"Loaded combined all-20 table: {combined.shape}")
    print(f"Loaded direct plasmid guide table: {plasmid.shape}")

    keep_cols = [
        "target_gene",
        "guide_id",
        "total_hit_rows",
        "direct_plasmid_contig_hits",
        "direct_chromosome_contig_hits",
        "no_fasta_header_match_hits",
        "unclassified_contig_hits",
        "direct_plasmid_hit_fraction",
        "direct_plasmid_mapping_decision",
    ]

    plasmid = plasmid[keep_cols].copy()

    merged = combined.merge(
        plasmid,
        on=["target_gene", "guide_id"],
        how="left"
    )

    zero_cols = [
        "total_hit_rows",
        "direct_plasmid_contig_hits",
        "direct_chromosome_contig_hits",
        "no_fasta_header_match_hits",
        "unclassified_contig_hits",
        "direct_plasmid_hit_fraction",
    ]

    for col in zero_cols:
        merged[col] = merged[col].fillna(0)

    int_cols = [
        "total_hit_rows",
        "direct_plasmid_contig_hits",
        "direct_chromosome_contig_hits",
        "no_fasta_header_match_hits",
        "unclassified_contig_hits",
    ]

    for col in int_cols:
        merged[col] = merged[col].astype(int)

    merged["direct_plasmid_mapping_decision"] = merged[
        "direct_plasmid_mapping_decision"
    ].fillna("Zero_hit_no_direct_plasmid_contig_evidence")

    def direct_plasmid_confidence(row):
        hits = row["direct_plasmid_contig_hits"]
        fraction = row["direct_plasmid_hit_fraction"]

        if hits == 0:
            return "No_direct_plasmid_evidence"
        if fraction >= 0.50:
            return "Strong_direct_plasmid_evidence"
        if hits >= 2:
            return "Moderate_direct_plasmid_evidence"
        return "Low_direct_plasmid_evidence_single_hit"

    merged["direct_plasmid_evidence_strength"] = merged.apply(
        direct_plasmid_confidence,
        axis=1
    )

    merged.to_csv(OUT_ALL20, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Direct Plasmid Evidence All-20 Summary\n")
        f.write("===================================================\n\n")
        f.write(f"Input combined all-20 table: {COMBINED_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Input direct plasmid guide table: {DIRECT_PLASMID_GUIDE.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output all-20 table: {OUT_ALL20.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Final guide rows: {len(merged)}\n")
        f.write(f"Unique guides: {merged['guide_id'].nunique()}\n")
        f.write(f"Zero-hit guides retained: {int((merged['total_hit_rows'] == 0).sum())}\n\n")

        f.write("Direct plasmid evidence strength counts:\n")
        for k, v in merged["direct_plasmid_evidence_strength"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nGuides with direct plasmid-contig evidence:\n")
        subset = merged[merged["direct_plasmid_contig_hits"] > 0]
        if subset.empty:
            f.write("- None\n")
        else:
            for _, r in subset.iterrows():
                f.write(
                    f"- {r['target_gene']} | {r['guide_id']} | "
                    f"plasmid_hits={int(r['direct_plasmid_contig_hits'])} | "
                    f"total_hits={int(r['total_hit_rows'])} | "
                    f"fraction={r['direct_plasmid_hit_fraction']:.3f} | "
                    f"strength={r['direct_plasmid_evidence_strength']}\n"
                )

    print(f"Wrote: {OUT_ALL20.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nDirect plasmid evidence strength counts:")
    print(merged["direct_plasmid_evidence_strength"].value_counts().to_string())

    print("\nRequirement 3 direct plasmid all-20 merge completed.")


if __name__ == "__main__":
    main()
