#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

COMBINED_TABLE = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_combined/requirement3_combined_requirement2_requirement3_final_guide_table_ALL20.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_combined"
OUT_REPORT = OUT_DIR / "requirement3_progress_report_text_for_docx.txt"


def main():
    if not COMBINED_TABLE.exists():
        raise FileNotFoundError(f"Combined table not found: {COMBINED_TABLE}")

    df = pd.read_csv(COMBINED_TABLE)

    top = (
        df.sort_values(["target_gene", "combined_req2_req3_rank_within_gene"])
        .groupby("target_gene")
        .head(2)
        .copy()
    )

    decision_counts = df["combined_req2_req3_final_decision"].value_counts()

    with open(OUT_REPORT, "w") as f:
        f.write("Additional Update: Requirement 3 Staged Mobile-Context Mapping Completed\n")
        f.write("=======================================================================\n\n")

        f.write("Purpose of this update\n")
        f.write("----------------------\n")
        f.write(
            "This update documents the staged completion of Requirement 3 for the "
            "FOTR-CRISPR computational manuscript: plasmid, resistance-island, SCCmec, "
            "integron, operon-proxy, and mobile-element context mapping. Requirement 1 "
            "had already produced the expanded Cas-OFFinder hit table, and Requirement 2 "
            "had already assigned database-refined functional annotations. Requirement 3 "
            "adds local genomic-neighborhood interpretation around each off-target hit.\n\n"
        )

        f.write("What was done\n")
        f.write("-------------\n")
        f.write(
            "Each expanded off-target hit was mapped to a local genomic neighborhood using "
            "the matched GFF/GFF3 annotations from the 36-accession expanded panel. For "
            "each hit, nearby annotated features within a +/-10 kb window were scanned for "
            "mobile-element, plasmid-associated, integron-like, SCCmec-like, AMR/resistance, "
            "and resistance-island-like signals. A simple operon-proxy flag was also generated "
            "based on nearby same-strand CDS features within a compact local region.\n\n"
        )

        f.write(
            "The hit-level context signals were then aggregated into guide-level mobile-context "
            "burden scores. Because hit-level aggregation naturally drops zero-hit guides, the "
            "Requirement 3 guide-level output was merged back into the all-20 Requirement 2 guide "
            "table. This restored the zero-hit guide and preserved all selected guides in the final "
            "audit table.\n\n"
        )

        f.write(
            "Finally, the Requirement 3 mobile-context burden was combined with the Requirement 2 "
            "database-refined functional burden. This created a combined final decision table that "
            "uses both biological functional risk and mobile-genetic-context risk. This prevents "
            "a guide from being incorrectly promoted only because it has low mobile-context burden "
            "when its broader functional off-target burden is high.\n\n"
        )

        f.write("Where the work was done\n")
        f.write("-----------------------\n")
        f.write("Project root:\n")
        f.write("- /Users/sanghati/Documents/crispr-mdr-resensitization\n\n")

        f.write("Main scripts:\n")
        f.write("- analysis/cas_offinder/requirement3_mobile_context_mapping.py\n")
        f.write("- analysis/cas_offinder/requirement3_merge_all20_guides.py\n")
        f.write("- analysis/cas_offinder/requirement3_make_combined_final_decision.py\n\n")

        f.write("Main output folder:\n")
        f.write("- results_cas_offinder/expanded_panel/mobile_context_mapping/\n\n")

        f.write("Final combined output folder:\n")
        f.write("- results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_combined/\n\n")

        f.write("Results obtained\n")
        f.write("----------------\n")
        f.write(f"Total final guide rows: {len(df)}\n")
        f.write(f"Unique selected guides retained: {df['guide_id'].nunique()}\n")
        f.write(f"Target genes represented: {df['target_gene'].nunique()}\n\n")

        f.write("Combined final decision counts:\n")
        for k, v in decision_counts.items():
            f.write(f"- {k}: {v}\n")

        f.write("\nFinal top 2 guides per gene after combined Requirement 2 + Requirement 3 scoring:\n")
        for _, r in top.iterrows():
            f.write(
                f"- {r['target_gene']} rank {int(r['combined_req2_req3_rank_within_gene'])}: "
                f"{r['guide_id']} | "
                f"Cas-OFFinder hits={int(r['total_hits_req3'])} | "
                f"Requirement 2 burden={r['db_refined_total_burden']} | "
                f"Requirement 3 mobile-context burden={r['requirement3_mobile_context_burden']} | "
                f"Final decision={r['combined_req2_req3_final_decision']}\n"
            )

        f.write("\nScientific interpretation\n")
        f.write("-------------------------\n")
        f.write(
            "Requirement 3 strengthens the FOTR-CRISPR framework by adding genomic-context "
            "interpretation to predicted off-target hits. Instead of only asking whether an "
            "off-target overlaps a coding or functional feature, this layer asks whether the "
            "off-target sits near mobile genetic elements, plasmid-associated genes, integron-like "
            "features, SCCmec-like regions, or AMR-associated neighborhoods. This is important for "
            "AMR-focused CRISPR guide design because many clinically relevant resistance genes are "
            "carried on plasmids, resistance islands, integrons, transposons, or SCCmec-like mobile "
            "genetic regions.\n\n"
        )

        f.write(
            "The final combined table shows that blaKPC_riskaware_top5 and blaNDM1_riskaware_top1 "
            "remain the cleanest computational candidates after both database-refined functional "
            "annotation and mobile-context mapping. The mcr1 and mecA candidates require more caution "
            "because their best remaining guides still show either elevated functional burden or "
            "mobile-context burden. This result is computational and should be described as staged "
            "candidate prioritization, not experimental validation.\n\n"
        )

        f.write("Safe manuscript wording\n")
        f.write("-----------------------\n")
        f.write(
            "We extended the expanded FOTR-CRISPR off-target annotation with local genomic-context "
            "mapping. For each candidate off-target hit, matched GFF/GFF3 annotations were scanned "
            "within a +/-10 kb neighborhood to identify mobile-element, plasmid-associated, "
            "integron-like, SCCmec-like, resistance-island-like, and operon-proxy signals. These "
            "signals were aggregated into guide-level mobile-context burden scores and merged back "
            "with the all-20 guide table to preserve zero-hit guides. The Requirement 3 context layer "
            "was then combined with Requirement 2 database-refined functional burden to generate a "
            "final guide-prioritization table. This staged analysis identified blaKPC_riskaware_top5 "
            "and blaNDM1_riskaware_top1 as the cleanest computational candidates, while mcr1 and mecA "
            "guides required caution due to elevated functional or mobile-context burden. These findings "
            "remain computational predictions and do not represent direct experimental validation or "
            "specialized-tool confirmation of plasmid type, SCCmec type, or integron structure.\n\n"
        )

        f.write("Current completion status\n")
        f.write("-------------------------\n")
        f.write(
            "Requirement 3 is now substantially completed as a staged computational layer. The current "
            "completion is approximately 70% because the pipeline now includes GFF-neighborhood-based "
            "mobile-context mapping, all-guide restoration, zero-hit retention, guide-level burden scoring, "
            "and combined Requirement 2 + Requirement 3 prioritization. The remaining optional extension "
            "would be direct specialized-tool confirmation using PlasmidFinder, SCCmecFinder, IntegronFinder, "
            "ISfinder/ISEScan, or MobileElementFinder.\n"
        )

    print(f"Wrote report text: {OUT_REPORT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
