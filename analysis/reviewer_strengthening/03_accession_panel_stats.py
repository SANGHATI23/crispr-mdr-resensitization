
import pandas as pd
import re
from pathlib import Path

OUT_DIR = Path("results_reviewer_strengthening/accession_panel_stats")
OUT_DIR.mkdir(parents=True, exist_ok=True)

HIT_FILE = Path("results_cas_offinder/expanded_panel/final_requirement1/requirement1_expanded_casoffinder_genomewide_offtarget_hit_table.csv")
GUIDE_FILE = Path("results_cas_offinder/expanded_panel/final_requirement1/requirement1_expanded_casoffinder_guide_level_summary.csv")
CONSERVATION_FILE = Path("results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/final_requirement4_complete/requirement4_final_integrated_conserved_offtarget_table_ALL20.csv")

def infer_species_or_strain(text):
    text = str(text)
    # Try to capture common phrase after accession in contig description
    # Example: Klebsiella pneumoniae strain ...
    m = re.search(r'(Klebsiella pneumoniae|Staphylococcus aureus|Escherichia coli|Acinetobacter baumannii|Pseudomonas aeruginosa|Enterococcus faecium|Enterococcus faecalis)', text, flags=re.I)
    if m:
        return m.group(1)
    return "Unclassified_from_header"

def infer_contig_class(text):
    t = str(text).lower()
    if "plasmid" in t:
        return "plasmid_labeled"
    if "chromosome" in t:
        return "chromosome_labeled"
    return "unclassified_contig_label"

def main():
    hits = pd.read_csv(HIT_FILE)
    guides = pd.read_csv(GUIDE_FILE)
    conservation = pd.read_csv(CONSERVATION_FILE)

    # Accessions represented in hit table
    accession_panel = (
        hits[["accession", "contig_description"]]
        .drop_duplicates()
        .copy()
    )

    accession_panel["inferred_species_from_header"] = accession_panel["contig_description"].apply(infer_species_or_strain)
    accession_panel["contig_label_class"] = accession_panel["contig_description"].apply(infer_contig_class)

    accession_summary = (
        accession_panel.groupby("accession")
        .agg(
            unique_contigs_with_hits=("contig_description", "nunique"),
            inferred_species_from_header=("inferred_species_from_header", lambda x: "; ".join(sorted(set(x)))),
            plasmid_labeled_contigs=("contig_label_class", lambda x: (x == "plasmid_labeled").sum()),
            chromosome_labeled_contigs=("contig_label_class", lambda x: (x == "chromosome_labeled").sum()),
            unclassified_contig_labels=("contig_label_class", lambda x: (x == "unclassified_contig_label").sum()),
        )
        .reset_index()
    )

    gene_accession = (
        hits.groupby(["target_gene", "accession"])
        .agg(
            hit_rows=("guide_id", "count"),
            guides_with_hits=("guide_id", "nunique"),
            unique_offtarget_sequences=("off_target_sequence", "nunique"),
            min_mismatch=("mismatch_count", "min"),
            max_mismatch=("mismatch_count", "max"),
        )
        .reset_index()
    )

    gene_panel_summary = (
        gene_accession.groupby("target_gene")
        .agg(
            accessions_with_any_hit=("accession", "nunique"),
            total_hit_rows=("hit_rows", "sum"),
            mean_hits_per_accession_with_hits=("hit_rows", "mean"),
            median_hits_per_accession_with_hits=("hit_rows", "median"),
            max_hits_in_single_accession=("hit_rows", "max"),
            total_guides_with_hits_sum=("guides_with_hits", "sum"),
        )
        .reset_index()
    )

    total_accessions_screened = int(conservation["total_accessions_screened"].max()) if "total_accessions_screened" in conservation.columns else accession_summary["accession"].nunique()

    overall = pd.DataFrame([{
        "total_accessions_screened_reported": total_accessions_screened,
        "accessions_with_any_offtarget_hit": accession_summary["accession"].nunique(),
        "unique_target_genes": hits["target_gene"].nunique(),
        "unique_guides_with_hits": hits["guide_id"].nunique(),
        "total_hit_rows": len(hits),
        "unique_contig_descriptions_with_hits": hits["contig_description"].nunique(),
        "plasmid_labeled_hit_contigs": int((accession_panel["contig_label_class"] == "plasmid_labeled").sum()),
        "chromosome_labeled_hit_contigs": int((accession_panel["contig_label_class"] == "chromosome_labeled").sum()),
        "unclassified_hit_contigs": int((accession_panel["contig_label_class"] == "unclassified_contig_label").sum()),
    }])

    accession_panel.to_csv(OUT_DIR / "accession_contig_panel_annotated.csv", index=False)
    accession_summary.to_csv(OUT_DIR / "accession_panel_summary.csv", index=False)
    gene_accession.to_csv(OUT_DIR / "gene_accession_hit_distribution.csv", index=False)
    gene_panel_summary.to_csv(OUT_DIR / "gene_panel_summary.csv", index=False)
    overall.to_csv(OUT_DIR / "overall_accession_panel_stats.csv", index=False)

    txt = []
    txt.append("Accession-panel statistics for expanded off-target screen")
    txt.append("========================================================")
    txt.append("")
    txt.append(f"Input hit table: {HIT_FILE}")
    txt.append(f"Reported total accessions screened: {total_accessions_screened}")
    txt.append(f"Accessions with at least one predicted off-target hit: {accession_summary['accession'].nunique()}")
    txt.append(f"Total parsed off-target hit rows: {len(hits)}")
    txt.append(f"Unique guides with at least one hit: {hits['guide_id'].nunique()}")
    txt.append(f"Unique target genes represented: {hits['target_gene'].nunique()}")
    txt.append("")
    txt.append("Overall panel stats:")
    txt.append(overall.to_string(index=False))
    txt.append("")
    txt.append("Gene-level accession distribution:")
    txt.append(gene_panel_summary.to_string(index=False))
    txt.append("")
    txt.append("Reviewer-facing interpretation:")
    txt.append("This analysis makes the 36-accession screening panel auditable by summarizing which accessions, contig labels, and target-gene hit distributions contributed to the expanded off-target screen.")
    txt.append("It does not claim exhaustive coverage of all clinical strain diversity. Instead, it supports a transparent computational validation panel and identifies how off-target burden is distributed across the screened genomes.")
    txt.append("")
    txt.append("Recommended manuscript wording:")
    txt.append("To make the expanded genome panel auditable, we summarized accession-level and gene-level hit distributions across the 36 screened genomes. The panel produced 1,997 parsed candidate off-target hits across all four target genes. Accession-level summaries were used to distinguish broad recurrent off-target burden from isolated accession-specific signals. This panel is not intended to exhaustively represent all clinical strain diversity, but it provides a transparent computational stress panel for comparing guide-level off-target recurrence.")

    (OUT_DIR / "accession_panel_stats_summary.txt").write_text("\n".join(txt))

    print("Wrote accession-panel stats outputs.")
    print()
    print("Overall:")
    print(overall.to_string(index=False))
    print()
    print("Gene-level summary:")
    print(gene_panel_summary.to_string(index=False))

if __name__ == "__main__":
    main()
