
import pandas as pd
from pathlib import Path

OUT_DIR = Path("manuscript_reviewer_strengthening")
OUT_DIR.mkdir(parents=True, exist_ok=True)

RGI_SUMMARY = Path("results_reviewer_strengthening/rgi_card/card_rgi_dictionary_confirmation_summary_by_gene.csv")
COMPARATOR_SUMMARY = Path("results_reviewer_strengthening/comparator_overlap/comparator_overlap_summary_cfd_rs3.csv")
CI_SUMMARY = Path("results_reviewer_strengthening/variance_ci/conserved_burden_variance_ci_by_guide.csv")
ACCESSION_OVERALL = Path("results_reviewer_strengthening/accession_panel_stats/overall_accession_panel_stats.csv")
ACCESSION_GENE = Path("results_reviewer_strengthening/accession_panel_stats/gene_panel_summary.csv")

def read_checked(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)

def main():
    rgi = read_checked(RGI_SUMMARY)
    comparator = read_checked(COMPARATOR_SUMMARY)
    ci = read_checked(CI_SUMMARY)
    accession_overall = read_checked(ACCESSION_OVERALL)
    accession_gene = read_checked(ACCESSION_GENE)

    # 1. Reviewer strengthening summary table
    strengthening_rows = [
        {
            "Reviewer concern addressed": "Functional off-target annotation may rely too heavily on GFF/GFF3 dictionary matching.",
            "Additional analysis added": "CARD/RGI-style direct-confirmation candidate set prepared from dictionary-flagged AMR, virulence, essential-like, mobile, and high-confidence functional hits.",
            "Main result": f"{int(rgi['total_hits'].sum()) if 'total_hits' in rgi.columns else 'NA'} total hits summarized by gene; direct-confirmation candidate files generated for high-interest functional hits.",
            "How to frame in manuscript": "This does not replace direct CARD/RGI execution, but makes the database-confirmation subset explicit and auditable."
        },
        {
            "Reviewer concern addressed": "FOTR ranking may be incremental without comparison to independent guide-scoring methods.",
            "Additional analysis added": "Orthogonal comparator overlap analysis using CFD specificity and RS3 activity rankings.",
            "Main result": "; ".join([
                f"{row['comparison']}: overlap {int(row['overlap_count'])}/{int(row['comparator_set_size'])}, Jaccard {row['jaccard']}"
                for _, row in comparator.iterrows()
            ]),
            "How to frame in manuscript": "FOTR shows limited/partial overlap with single-objective comparators, supporting that it is a distinct multi-objective prioritization framework rather than a replica of CFD or RS3."
        },
        {
            "Reviewer concern addressed": "Final conserved off-target burden is reported as point estimates only.",
            "Additional analysis added": "Per-guide accession-level variance and approximate 95% confidence intervals for database-refined burden.",
            "Main result": "Clean candidates showed low accession-level burden, e.g., blaNDM1_riskaware_top1 mean 0.0833 with CI 0–0.2467; high-risk guides such as mecA_riskaware_top1 and mcr1_riskaware_top5 showed much wider/high burden intervals.",
            "How to frame in manuscript": "This checks whether guide burden is consistently low across accessions or concentrated in recurrent/high-burden accession-level signals."
        },
        {
            "Reviewer concern addressed": "The 36-accession genome panel may need clearer characterization.",
            "Additional analysis added": "Accession-panel statistics summarizing accessions, contigs, plasmid/chromosome labels, and gene-level hit distribution.",
            "Main result": f"{int(accession_overall.loc[0, 'total_accessions_screened_reported'])} accessions screened; {int(accession_overall.loc[0, 'total_hit_rows'])} hit rows; {int(accession_overall.loc[0, 'plasmid_labeled_hit_contigs'])} plasmid-labeled hit contigs; {int(accession_overall.loc[0, 'chromosome_labeled_hit_contigs'])} chromosome-labeled hit contigs.",
            "How to frame in manuscript": "The panel is not exhaustive clinical diversity, but it is now transparent and auditable as a computational stress panel."
        },
    ]

    summary_table = pd.DataFrame(strengthening_rows)
    summary_table.to_csv(OUT_DIR / "reviewer_strengthening_summary_table.csv", index=False)

    # 2. Manuscript-ready concise text
    text = []
    text.append("Reviewer-strengthening additions for FOTR-CRISPR")
    text.append("================================================")
    text.append("")
    text.append("Four additional computational checks were added to strengthen the manuscript against likely reviewer concerns.")
    text.append("")
    text.append("First, because the functional off-target layer originally used GFF/GFF3-derived annotations with curated dictionary refinement, we generated a CARD/RGI-style direct-confirmation candidate set. This table identifies the subset of off-target hits most suitable for database-level confirmation, including AMR/resistance-related, virulence-related, essential-like, mobile-element-related, and high-confidence functional hits.")
    text.append("")
    text.append("Second, to test whether FOTR simply reproduces existing single-objective scoring behavior, we compared the final FOTR-selected guides against orthogonal CFD specificity and RS3 activity rankings. The comparator overlap was limited rather than complete: FOTR top-20 versus CFD top-20 overlapped by 9 guides with Jaccard 0.2903, while FOTR top-20 versus RS3 top-20 showed no overlap in the available comparison table. This should be framed carefully as evidence that FOTR is a distinct multi-objective prioritization framework, not as strong agreement with external tools.")
    text.append("")
    text.append("Third, to avoid reporting conserved off-target burden as point estimates only, we computed accession-level burden variance and approximate 95% confidence intervals for each guide. Clean candidates such as blaNDM1_riskaware_top1 showed very low mean accession burden and narrow confidence intervals, whereas high-concern guides such as mecA_riskaware_top1 and mcr1_riskaware_top5 showed substantially higher and wider burden intervals.")
    text.append("")
    text.append("Fourth, to make the 36-accession screen more transparent, we added accession-panel statistics. The expanded screen covered 36 reported accessions and produced 1,997 parsed candidate off-target hit rows across four target genes. The hit-containing contig set included 13 plasmid-labeled and 33 chromosome-labeled contigs, with mcr1 and mecA showing the highest gene-level hit burden.")
    text.append("")
    text.append("Suggested manuscript insertion:")
    text.append("")
    text.append("To strengthen the computational validation layer, we added four reviewer-oriented checks. First, high-interest dictionary-flagged off-target hits were exported as a CARD/RGI-style direct-confirmation candidate set. Second, FOTR-selected guides were compared against independent CFD specificity and RS3 activity rankings. Limited overlap with these single-objective comparators supports that FOTR does not simply reproduce one existing scoring dimension, but instead integrates activity, specificity, conservation, and functional safety evidence. Third, accession-level variance and approximate 95% confidence intervals were computed for conserved off-target burden, distinguishing consistently low-burden candidates from recurrent high-burden guides. Fourth, accession-panel statistics were generated to make the 36-genome screen auditable. Together, these additions strengthen the manuscript by converting likely reviewer concerns into explicit computational checks while preserving the claim boundary that FOTR is a computational prioritization framework, not experimental validation.")

    (OUT_DIR / "reviewer_strengthening_manuscript_text.txt").write_text("\n".join(text))

    # 3. Compact table for paper
    compact_rows = [
        {
            "Added check": "CARD/RGI-style confirmation subset",
            "Purpose": "Identify dictionary-flagged functional hits suitable for direct database confirmation",
            "Key output": "High-interest candidate hit table by guide/gene",
            "Reviewer value": "Makes annotation limitation auditable"
        },
        {
            "Added check": "CFD/RS3 comparator overlap",
            "Purpose": "Compare FOTR with orthogonal specificity/activity rankings",
            "Key output": "CFD top-20 overlap 9/20, Jaccard 0.2903; RS3 top-20 overlap 0/20",
            "Reviewer value": "Shows FOTR is not just a single-objective score replica"
        },
        {
            "Added check": "Burden variance/CI",
            "Purpose": "Quantify accession-level stability of conserved burden",
            "Key output": "Mean, SD, SE, and 95% CI per guide",
            "Reviewer value": "Adds uncertainty/stability interpretation"
        },
        {
            "Added check": "Accession-panel statistics",
            "Purpose": "Characterize the 36-genome off-target screen",
            "Key output": "1,997 hits; 49 contig descriptions; 13 plasmid-labeled and 33 chromosome-labeled hit contigs",
            "Reviewer value": "Clarifies panel scope and burden distribution"
        },
    ]

    compact = pd.DataFrame(compact_rows)
    compact.to_csv(OUT_DIR / "manuscript_ready_reviewer_strengthening_compact_table.csv", index=False)

    print("Wrote:")
    print(OUT_DIR / "reviewer_strengthening_summary_table.csv")
    print(OUT_DIR / "reviewer_strengthening_manuscript_text.txt")
    print(OUT_DIR / "manuscript_ready_reviewer_strengthening_compact_table.csv")
    print()
    print(compact.to_string(index=False))

if __name__ == "__main__":
    main()
