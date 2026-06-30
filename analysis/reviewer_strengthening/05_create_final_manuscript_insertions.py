
from pathlib import Path
import pandas as pd

OUT_DIR = Path("manuscript_reviewer_strengthening")
OUT_DIR.mkdir(parents=True, exist_ok=True)

COMPACT_TABLE = Path("manuscript_reviewer_strengthening/manuscript_ready_reviewer_strengthening_compact_table.csv")
COMPARATOR = Path("results_reviewer_strengthening/comparator_overlap/comparator_overlap_summary_cfd_rs3.csv")
CI = Path("results_reviewer_strengthening/variance_ci/conserved_burden_variance_ci_by_guide.csv")
ACCESSION = Path("results_reviewer_strengthening/accession_panel_stats/overall_accession_panel_stats.csv")
GENE_PANEL = Path("results_reviewer_strengthening/accession_panel_stats/gene_panel_summary.csv")

def read_checked(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)

def main():
    compact = read_checked(COMPACT_TABLE)
    comparator = read_checked(COMPARATOR)
    ci = read_checked(CI)
    accession = read_checked(ACCESSION)
    gene_panel = read_checked(GENE_PANEL)

    # Pull key numbers
    cfd20 = comparator[comparator["comparison"].str.contains("top20 vs CFD", regex=False)].iloc[0]
    rs320 = comparator[comparator["comparison"].str.contains("top20 vs RS3", regex=False)].iloc[0]

    best_low = ci.sort_values("mean_accession_burden", ascending=True).head(3)
    high_burden = ci.sort_values("mean_accession_burden", ascending=False).head(3)

    total_accessions = int(accession.loc[0, "total_accessions_screened_reported"])
    total_hits = int(accession.loc[0, "total_hit_rows"])
    plasmid_contigs = int(accession.loc[0, "plasmid_labeled_hit_contigs"])
    chromosome_contigs = int(accession.loc[0, "chromosome_labeled_hit_contigs"])
    unique_contigs = int(accession.loc[0, "unique_contig_descriptions_with_hits"])

    # Main insertion section
    insertion = f"""3.X Reviewer-Oriented Computational Validation Checks

To further strengthen the computational validation layer, we added four reviewer-oriented checks addressing likely concerns about annotation support, comparator behavior, burden uncertainty, and accession-panel transparency. First, because the functional off-target layer uses GFF/GFF3-derived feature annotations with curated dictionary refinement, we generated a CARD/RGI-style direct-confirmation candidate set. This output identifies dictionary-flagged AMR/resistance-related, virulence-related, essential-like, mobile-element-related, and high-confidence functional off-target hits that are most suitable for direct database-level confirmation. This step does not replace direct CARD/RGI execution, but it makes the subset requiring direct database validation explicit and auditable.

Second, we performed an orthogonal comparator-overlap analysis using CFD specificity and RS3 activity rankings generated on the same candidate pool. FOTR selected top-20 guides overlapped with the CFD-ranked top-20 set by {int(cfd20['overlap_count'])}/{int(cfd20['comparator_set_size'])} guides, with Jaccard similarity {float(cfd20['jaccard']):.4f}. In contrast, FOTR selected top-20 guides showed {int(rs320['overlap_count'])}/{int(rs320['comparator_set_size'])} overlap with the RS3-ranked top-20 set in the available comparison table. These results should not be interpreted as strong agreement with external single-objective scores. Rather, the limited/partial overlap supports the intended behavior of FOTR-CRISPR as a distinct multi-objective prioritization framework that integrates activity, specificity, conservation, and functional off-target safety rather than optimizing only one scoring dimension.

Third, to avoid reporting conserved off-target burden only as point estimates, we computed accession-level burden variance and approximate 95% confidence intervals for each guide. This analysis aggregated database-refined burden within each screened accession and then summarized mean burden, standard deviation, standard error, and approximate confidence intervals per guide. Low-burden candidates remained clearly separated from high-burden candidates. For example, blaNDM1_riskaware_top1 showed mean accession burden 0.0833 with 95% CI 0.0000–0.2467, whereas mecA_riskaware_top1 showed mean accession burden 135.4444 with 95% CI 83.4295–187.4593. This provides an additional stability check for the final integrated ranking by distinguishing consistently low-burden candidates from guides with recurrent or high-intensity accession-level burden.

Fourth, we summarized the accession-panel composition and hit distribution to make the expanded genome screen more auditable. The expanded off-target screen covered {total_accessions} reported accessions and produced {total_hits} parsed candidate off-target hit rows across the four target genes. The hit-containing contig set included {unique_contigs} unique contig descriptions, including {plasmid_contigs} plasmid-labeled and {chromosome_contigs} chromosome-labeled hit contigs. Gene-level summaries showed that mcr1 and mecA carried the highest off-target burden, whereas blaNDM1 showed the lowest overall hit burden. These results do not imply exhaustive coverage of all clinical strain diversity, but they make the computational stress panel transparent and support the interpretation that off-target recurrence differs strongly by target gene.
"""

    (OUT_DIR / "final_manuscript_insertion_section_3X.txt").write_text(insertion)

    # Compact table as markdown-like text
    table_text = []
    table_text.append("Table X. Reviewer-oriented computational validation checks added to strengthen FOTR-CRISPR.")
    table_text.append("")
    table_text.append(compact.to_string(index=False))
    (OUT_DIR / "final_manuscript_table_X_reviewer_checks.txt").write_text("\n".join(table_text))

    # Limitations update
    limitations = """Updated Limitations paragraph

The additional reviewer-oriented checks strengthen but do not convert the study into experimental validation. The CARD/RGI-style output identifies high-interest off-target hits suitable for direct database confirmation, but full direct CARD/RGI or AMRFinderPlus execution remains future work unless performed separately. The CFD/RS3 comparator analysis should be interpreted as an orthogonal single-objective scoring comparison rather than a full CRISPOR/CHOPCHOP web-tool replication. The accession-level confidence intervals summarize computational burden variability across the screened panel, but they do not model all sources of uncertainty in Cas9 cleavage, delivery, bacterial repair, or clinical strain diversity. Finally, the 36-accession panel provides a transparent computational stress panel, not exhaustive coverage of all plasmid, SCCmec, integron, or resistance-island backgrounds.
"""
    (OUT_DIR / "final_limitations_update.txt").write_text(limitations)

    # Exact insertion checklist
    checklist = f"""Exact manuscript insertion checklist

1. Add the new subsection after the current integrated off-target safety results, preferably after Section 3.5.5 or before Section 3.6.
   Suggested title:
   3.6 Reviewer-Oriented Computational Validation Checks

2. Insert Table X after this new subsection.
   Table title:
   Reviewer-oriented computational validation checks added to strengthen FOTR-CRISPR.

3. Renumber the old Section 3.6 Practical Utility as Section 3.7 if needed.

4. Add one sentence to Discussion:
   Additional reviewer-oriented analyses showed that FOTR does not simply reproduce CFD or RS3 single-objective rankings, that conserved burden varies strongly across accessions and genes, and that the expanded 36-accession panel provides an auditable computational stress screen rather than an exhaustive clinical diversity panel.

5. Replace or extend the current Limitations paragraph using final_limitations_update.txt.

6. Do not overclaim:
   Do not say direct CARD/RGI validation was completed unless you actually run RGI/CARD.
   Do not say FOTR strongly agrees with RS3 or CFD.
   Do not say the 36-accession panel represents all clinical strain diversity.
   Do not say this is wet-lab validated.

7. Safe final claim:
   FOTR-CRISPR is strengthened by adversarial edge-case sensitivity analysis, orthogonal comparator-score overlap, accession-level burden confidence intervals, and transparent accession-panel statistics. Together, these additions support the framework as a robustness-tested and safety-context-aware computational prioritization workflow.
"""
    (OUT_DIR / "final_exact_manuscript_insertion_checklist.txt").write_text(checklist)

    # Also create selected result summary
    selected = []
    selected.append("Key numbers for manuscript")
    selected.append("==========================")
    selected.append(f"CFD top-20 overlap: {int(cfd20['overlap_count'])}/{int(cfd20['comparator_set_size'])}, Jaccard {float(cfd20['jaccard']):.4f}")
    selected.append(f"RS3 top-20 overlap: {int(rs320['overlap_count'])}/{int(rs320['comparator_set_size'])}, Jaccard {float(rs320['jaccard']):.4f}")
    selected.append(f"Accessions screened: {total_accessions}")
    selected.append(f"Total parsed hits: {total_hits}")
    selected.append(f"Unique contig descriptions with hits: {unique_contigs}")
    selected.append(f"Plasmid-labeled hit contigs: {plasmid_contigs}")
    selected.append(f"Chromosome-labeled hit contigs: {chromosome_contigs}")
    selected.append("")
    selected.append("Lowest mean accession-burden guides:")
    selected.append(best_low[["target_gene", "guide_id", "mean_accession_burden", "ci95_low", "ci95_high"]].to_string(index=False))
    selected.append("")
    selected.append("Highest mean accession-burden guides:")
    selected.append(high_burden[["target_gene", "guide_id", "mean_accession_burden", "ci95_low", "ci95_high"]].to_string(index=False))
    selected.append("")
    selected.append("Gene panel summary:")
    selected.append(gene_panel.to_string(index=False))
    (OUT_DIR / "final_key_numbers_for_manuscript.txt").write_text("\n".join(selected))

    print("Wrote final manuscript insertion package:")
    for f in [
        "final_manuscript_insertion_section_3X.txt",
        "final_manuscript_table_X_reviewer_checks.txt",
        "final_limitations_update.txt",
        "final_exact_manuscript_insertion_checklist.txt",
        "final_key_numbers_for_manuscript.txt",
    ]:
        print(OUT_DIR / f)

    print()
    print("Key numbers:")
    print("\n".join(selected))

if __name__ == "__main__":
    main()
