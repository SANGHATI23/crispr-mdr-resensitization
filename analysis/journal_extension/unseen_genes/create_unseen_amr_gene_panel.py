from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]

OUT_DIR = ROOT / "results_journal_extension" / "unseen_genes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "unseen_amr_gene_panel.csv"
OUT_SUMMARY = OUT_DIR / "unseen_amr_gene_panel_summary.txt"

GENES = [
    {
        "gene": "vanA",
        "resistance_category": "glycopeptide_resistance",
        "primary_drug_class": "vancomycin",
        "common_context": "vancomycin-resistant Enterococcus and transferable van gene clusters",
        "mechanism_summary": "Alters peptidoglycan precursor termini to reduce vancomycin binding.",
        "extension_role": "Unseen AMR gene for prospective validation of FOTR-CRISPR generalizability.",
        "clinical_priority": "High",
        "expected_target_context": "mobile_element_or_gene_cluster",
    },
    {
        "gene": "tetM",
        "resistance_category": "tetracycline_resistance",
        "primary_drug_class": "tetracycline",
        "common_context": "conjugative transposons and mobile resistance elements",
        "mechanism_summary": "Ribosomal protection protein that reduces tetracycline inhibition.",
        "extension_role": "Unseen AMR gene with mobile-element relevance.",
        "clinical_priority": "Moderate",
        "expected_target_context": "mobile_element_or_conjugative_element",
    },
    {
        "gene": "ermB",
        "resistance_category": "macrolide_lincosamide_streptogramin_resistance",
        "primary_drug_class": "macrolide",
        "common_context": "Gram-positive pathogens and mobile resistance elements",
        "mechanism_summary": "rRNA methyltransferase that reduces macrolide binding to the ribosome.",
        "extension_role": "Unseen AMR gene for testing functional-context scoring outside beta-lactam/colistin targets.",
        "clinical_priority": "Moderate",
        "expected_target_context": "mobile_element_or_resistance_gene_cluster",
    },
    {
        "gene": "aac6Ib",
        "resistance_category": "aminoglycoside_resistance",
        "primary_drug_class": "aminoglycoside",
        "common_context": "plasmids, integrons, and multidrug resistance backgrounds",
        "mechanism_summary": "Aminoglycoside acetyltransferase that enzymatically modifies aminoglycosides.",
        "extension_role": "Unseen AMR gene for plasmid/integron-associated resistance testing.",
        "clinical_priority": "High",
        "expected_target_context": "plasmid_or_integron",
    },
    {
        "gene": "qnrS",
        "resistance_category": "fluoroquinolone_resistance",
        "primary_drug_class": "fluoroquinolone",
        "common_context": "plasmid-mediated quinolone resistance",
        "mechanism_summary": "Protects DNA gyrase/topoisomerase IV from quinolone inhibition.",
        "extension_role": "Unseen AMR gene for plasmid-mediated resistance validation.",
        "clinical_priority": "Moderate",
        "expected_target_context": "plasmid",
    },
]


def main():
    df = pd.DataFrame(GENES)
    df["panel_status"] = "prospective_unseen_gene"
    df["used_in_original_bibm_core"] = "No"

    df.to_csv(OUT_TABLE, index=False)

    lines = []
    lines.append("Prospective Unseen AMR Gene Panel Summary")
    lines.append("=========================================")
    lines.append("")
    lines.append(f"Output table: {OUT_TABLE}")
    lines.append(f"Total unseen genes: {len(df)}")
    lines.append("")
    lines.append("Genes included:")
    for _, r in df.iterrows():
        lines.append(
            f"- {r['gene']}: {r['primary_drug_class']} resistance | "
            f"priority={r['clinical_priority']} | context={r['expected_target_context']}"
        )
    lines.append("")
    lines.append("Interpretation:")
    lines.append(
        "This panel defines prospective AMR targets not used in the original four-gene FOTR-CRISPR core. "
        "It is intended for journal-extension generalizability testing before guide generation, external-genome screening, "
        "and multiplex guide-set design."
    )

    OUT_SUMMARY.write_text("\n".join(lines))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
