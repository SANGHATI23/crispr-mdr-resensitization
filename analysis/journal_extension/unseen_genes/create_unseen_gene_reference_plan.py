from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]

INPUT_PANEL = ROOT / "results_journal_extension" / "unseen_genes" / "unseen_amr_gene_panel.csv"

OUT_DIR = ROOT / "results_journal_extension" / "unseen_genes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "unseen_gene_reference_plan.csv"
OUT_SUMMARY = OUT_DIR / "unseen_gene_reference_plan_summary.txt"

REFERENCE_PLAN = {
    "vanA": {
        "preferred_reference_source": "CARD or NCBI RefSeq",
        "target_sequence_type": "coding_sequence",
        "recommended_query": "vanA vancomycin resistance gene complete coding sequence Enterococcus",
        "expected_host_context": "Enterococcus faecium / Enterococcus faecalis",
        "collection_priority": 1,
        "notes": "Use canonical vanA CDS first; later add variant sequences from VRE genomes."
    },
    "tetM": {
        "preferred_reference_source": "CARD or NCBI RefSeq",
        "target_sequence_type": "coding_sequence",
        "recommended_query": "tetM tetracycline resistance gene complete coding sequence",
        "expected_host_context": "Streptococcus / Enterococcus / mobile conjugative elements",
        "collection_priority": 2,
        "notes": "Use representative tetM CDS; later include mobile-element-associated variants."
    },
    "ermB": {
        "preferred_reference_source": "CARD or NCBI RefSeq",
        "target_sequence_type": "coding_sequence",
        "recommended_query": "ermB macrolide resistance gene complete coding sequence",
        "expected_host_context": "Streptococcus / Enterococcus / Staphylococcus",
        "collection_priority": 3,
        "notes": "Use representative ermB CDS for unseen-gene guide enumeration."
    },
    "aac6Ib": {
        "preferred_reference_source": "CARD or NCBI RefSeq",
        "target_sequence_type": "coding_sequence",
        "recommended_query": "aac(6')-Ib aminoglycoside acetyltransferase complete coding sequence",
        "expected_host_context": "Enterobacteriaceae plasmids and integrons",
        "collection_priority": 4,
        "notes": "Normalize gene name as aac6Ib for filenames; preserve biological label as aac(6')-Ib in manuscript."
    },
    "qnrS": {
        "preferred_reference_source": "CARD or NCBI RefSeq",
        "target_sequence_type": "coding_sequence",
        "recommended_query": "qnrS quinolone resistance gene complete coding sequence",
        "expected_host_context": "plasmid-mediated quinolone resistance Enterobacteriaceae",
        "collection_priority": 5,
        "notes": "Use qnrS representative CDS; later test plasmid-background conservation."
    },
}


def main():
    panel = pd.read_csv(INPUT_PANEL)

    rows = []
    for _, r in panel.iterrows():
        gene = r["gene"]
        p = REFERENCE_PLAN.get(gene, {})
        rows.append({
            "gene": gene,
            "resistance_category": r.get("resistance_category", ""),
            "primary_drug_class": r.get("primary_drug_class", ""),
            "clinical_priority": r.get("clinical_priority", ""),
            "preferred_reference_source": p.get("preferred_reference_source", "CARD or NCBI RefSeq"),
            "target_sequence_type": p.get("target_sequence_type", "coding_sequence"),
            "recommended_query": p.get("recommended_query", ""),
            "expected_host_context": p.get("expected_host_context", ""),
            "collection_priority": p.get("collection_priority", ""),
            "planned_fasta_path": f"data_journal_extension/unseen_genes/{gene}.fa",
            "planned_guides_output": f"results_journal_extension/unseen_genes/{gene}_candidate_guides.csv",
            "notes": p.get("notes", "")
        })

    df = pd.DataFrame(rows).sort_values("collection_priority")
    df.to_csv(OUT_TABLE, index=False)

    lines = []
    lines.append("Unseen Gene Reference Plan Summary")
    lines.append("==================================")
    lines.append("")
    lines.append(f"Input panel: {INPUT_PANEL}")
    lines.append(f"Output table: {OUT_TABLE}")
    lines.append(f"Total genes planned: {len(df)}")
    lines.append("")
    lines.append("Planned sequence collection:")
    for _, r in df.iterrows():
        lines.append(
            f"- {r['gene']}: source={r['preferred_reference_source']} | "
            f"type={r['target_sequence_type']} | planned_fasta={r['planned_fasta_path']}"
        )
    lines.append("")
    lines.append("Interpretation:")
    lines.append(
        "This table defines the sequence-collection plan for prospective unseen AMR genes. "
        "It is intentionally separated from the guide-generation step so that sequence provenance, "
        "gene naming, and target scope can be audited before adding FASTA inputs."
    )

    OUT_SUMMARY.write_text("\n".join(lines))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
