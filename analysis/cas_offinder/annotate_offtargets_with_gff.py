import pandas as pd
from pathlib import Path
from urllib.parse import unquote

PROJECT_ROOT = Path(__file__).resolve().parents[2]

HITS_FILE = PROJECT_ROOT / "results_cas_offinder" / "parsed" / "cas_offinder_all_genome_offtarget_hits_parsed.csv"

ANNOTATION_DIR = PROJECT_ROOT / "data_cas_offinder" / "annotations"

OUTPUT_DIR = PROJECT_ROOT / "results_cas_offinder" / "annotated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "cas_offinder_offtarget_hits_gff_annotated.csv"
SUMMARY_FILE = OUTPUT_DIR / "cas_offinder_gff_annotation_summary.txt"

GFF_FILES = {
    "Klebsiella_pneumoniae_reference": ANNOTATION_DIR / "Klebsiella_pneumoniae_reference.gff3",
    "Ecoli_reference": ANNOTATION_DIR / "Ecoli_reference.gff3",
    "Staphylococcus_aureus_MRSA_reference": ANNOTATION_DIR / "Staphylococcus_aureus_MRSA_reference.gff3",
}


def parse_attributes(attr_string):
    """
    Parse GFF3 attributes like:
    ID=gene-abc;Name=abc;gene=abc;locus_tag=XYZ;product=...
    """
    attrs = {}
    if pd.isna(attr_string):
        return attrs

    for item in str(attr_string).split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = unquote(value.strip())
    return attrs


def load_gff_features(gff_path, genome_label):
    rows = []

    if not gff_path.exists():
        print(f"WARNING: Missing GFF3 for {genome_label}: {gff_path}")
        return pd.DataFrame()

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Keep biologically useful feature types.
            # CDS is best for protein-coding overlap.
            # gene captures gene intervals even if CDS formatting is inconsistent.
            if feature_type not in ["gene", "CDS", "pseudogene", "ncRNA", "rRNA", "tRNA"]:
                continue

            attrs = parse_attributes(attributes)

            gene_name = (
                attrs.get("gene")
                or attrs.get("Name")
                or attrs.get("locus_tag")
                or attrs.get("ID")
                or ""
            )

            product = attrs.get("product", "")

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            rows.append(
                {
                    "genome_label": genome_label,
                    "seqid": seqid,
                    "feature_type": feature_type,
                    "feature_start": start_i,
                    "feature_end": end_i,
                    "feature_strand": strand,
                    "gene_name": gene_name,
                    "product": product,
                    "raw_attributes": attributes,
                }
            )

    df = pd.DataFrame(rows)
    print(f"Loaded {genome_label} GFF features: {df.shape}")
    return df


def infer_seqid_from_contig_description(contig_description):
    """
    Cas-OFFinder contig_description begins with accession, e.g.
    NC_016845.1 Klebsiella pneumoniae...
    """
    if pd.isna(contig_description):
        return ""
    return str(contig_description).split()[0]


def annotate_hit(hit, feature_df):
    pos = hit["position"]
    seqid = hit["seqid"]

    if pd.isna(pos) or not seqid:
        return {
            "genomic_context": "Unknown",
            "overlapping_feature_type": "",
            "overlapping_gene_name": "",
            "overlapping_product": "",
            "overlapping_feature_start": None,
            "overlapping_feature_end": None,
            "overlapping_feature_strand": "",
        }

    overlapping = feature_df[
        (feature_df["seqid"] == seqid)
        & (feature_df["feature_start"] <= pos)
        & (feature_df["feature_end"] >= pos)
    ].copy()

    if overlapping.empty:
        return {
            "genomic_context": "Intergenic_or_unannotated",
            "overlapping_feature_type": "",
            "overlapping_gene_name": "",
            "overlapping_product": "",
            "overlapping_feature_start": None,
            "overlapping_feature_end": None,
            "overlapping_feature_strand": "",
        }

    # Prefer CDS over gene if both overlap.
    priority = {
        "CDS": 1,
        "gene": 2,
        "pseudogene": 3,
        "ncRNA": 4,
        "rRNA": 5,
        "tRNA": 6,
    }
    overlapping["priority"] = overlapping["feature_type"].map(priority).fillna(99)
    best = overlapping.sort_values(["priority", "feature_start"]).iloc[0]

    if best["feature_type"] == "CDS":
        context = "Coding_CDS"
    elif best["feature_type"] == "gene":
        context = "Gene_region_nonCDS_or_gene_feature"
    elif best["feature_type"] in ["ncRNA", "rRNA", "tRNA"]:
        context = "Noncoding_RNA"
    elif best["feature_type"] == "pseudogene":
        context = "Pseudogene"
    else:
        context = "Other_annotated_feature"

    return {
        "genomic_context": context,
        "overlapping_feature_type": best["feature_type"],
        "overlapping_gene_name": best["gene_name"],
        "overlapping_product": best["product"],
        "overlapping_feature_start": best["feature_start"],
        "overlapping_feature_end": best["feature_end"],
        "overlapping_feature_strand": best["feature_strand"],
    }


def main():
    print("Annotating Cas-OFFinder off-target hits with GFF3 features...")

    if not HITS_FILE.exists():
        raise FileNotFoundError(f"Missing parsed hits file: {HITS_FILE}")

    hits = pd.read_csv(HITS_FILE)
    print(f"Loaded parsed hits: {hits.shape}")

    hits["seqid"] = hits["contig_description"].apply(infer_seqid_from_contig_description)

    all_feature_tables = {}
    for genome_label, gff_path in GFF_FILES.items():
        all_feature_tables[genome_label] = load_gff_features(gff_path, genome_label)

    annotated_rows = []

    for _, hit in hits.iterrows():
        genome_label = hit["genome_label"]
        feature_df = all_feature_tables.get(genome_label, pd.DataFrame())

        if feature_df.empty:
            ann = {
                "genomic_context": "No_GFF_available",
                "overlapping_feature_type": "",
                "overlapping_gene_name": "",
                "overlapping_product": "",
                "overlapping_feature_start": None,
                "overlapping_feature_end": None,
                "overlapping_feature_strand": "",
            }
        else:
            ann = annotate_hit(hit, feature_df)

        row = hit.to_dict()
        row.update(ann)
        annotated_rows.append(row)

    annotated = pd.DataFrame(annotated_rows)

    # Simple first-pass biological flags from product/gene text.
    text = (
        annotated["overlapping_gene_name"].fillna("").astype(str)
        + " "
        + annotated["overlapping_product"].fillna("").astype(str)
    ).str.lower()

    annotated["keyword_amr_like"] = text.str.contains(
        "beta-lactamase|lactamase|resistance|efflux|multidrug|drug|antibiotic|methicillin|carbapenem|colistin|transposase|integrase",
        regex=True,
    )

    annotated["keyword_mobile_element_like"] = text.str.contains(
        "transposase|integrase|recombinase|insertion sequence|is[0-9]|phage|plasmid|mobile",
        regex=True,
    )

    annotated["keyword_essential_like"] = text.str.contains(
        "ribosomal|dna polymerase|rna polymerase|gyrase|topoisomerase|cell division|fts|replication|translation|elongation factor",
        regex=True,
    )

    annotated["first_pass_functional_offtarget_class"] = "Low_or_unknown_context"

    annotated.loc[
        annotated["genomic_context"].isin(["Intergenic_or_unannotated", "Unknown", "No_GFF_available"]),
        "first_pass_functional_offtarget_class",
    ] = "Intergenic_or_unannotated"

    annotated.loc[
        annotated["genomic_context"].isin(["Coding_CDS", "Gene_region_nonCDS_or_gene_feature"]),
        "first_pass_functional_offtarget_class",
    ] = "Coding_or_gene_overlap"

    annotated.loc[
        annotated["keyword_amr_like"],
        "first_pass_functional_offtarget_class",
    ] = "AMR_or_resistance_related_keyword"

    annotated.loc[
        annotated["keyword_mobile_element_like"],
        "first_pass_functional_offtarget_class",
    ] = "Mobile_element_related_keyword"

    annotated.loc[
        annotated["keyword_essential_like"],
        "first_pass_functional_offtarget_class",
    ] = "Essential_like_keyword"

    annotated.to_csv(OUTPUT_FILE, index=False)

    # Summaries
    context_summary = annotated.groupby("genomic_context").size().reset_index(name="hit_count")
    gene_summary = annotated.groupby(["mlcb_gene", "genomic_context"]).size().reset_index(name="hit_count")
    class_summary = annotated.groupby("first_pass_functional_offtarget_class").size().reset_index(name="hit_count")
    mismatch_context = (
        annotated.groupby(["mismatch_count", "genomic_context"])
        .size()
        .reset_index(name="hit_count")
        .sort_values(["mismatch_count", "genomic_context"])
    )

    high_interest = annotated[
        (annotated["mismatch_count"] <= 3)
        | (annotated["first_pass_functional_offtarget_class"].isin(
            [
                "AMR_or_resistance_related_keyword",
                "Mobile_element_related_keyword",
                "Essential_like_keyword",
            ]
        ))
    ].copy()

    high_interest_file = OUTPUT_DIR / "cas_offinder_high_interest_annotated_hits.csv"
    high_interest.to_csv(high_interest_file, index=False)

    with open(SUMMARY_FILE, "w") as f:
        f.write("GFF3 Functional Annotation Summary\n")
        f.write("=================================\n\n")
        f.write(f"Input hits: {HITS_FILE}\n")
        f.write(f"Annotated output: {OUTPUT_FILE}\n")
        f.write(f"Total hits annotated: {annotated.shape[0]}\n")
        f.write(f"High-interest hits: {high_interest.shape[0]}\n\n")

        f.write("Genomic context summary:\n")
        f.write(context_summary.to_string(index=False))
        f.write("\n\n")

        f.write("First-pass functional off-target class summary:\n")
        f.write(class_summary.to_string(index=False))
        f.write("\n\n")

        f.write("Gene x genomic context summary:\n")
        f.write(gene_summary.to_string(index=False))
        f.write("\n\n")

        f.write("Mismatch x genomic context summary:\n")
        f.write(mismatch_context.to_string(index=False))
        f.write("\n\n")

    print("\nWrote:")
    print(f"- {OUTPUT_FILE}")
    print(f"- {high_interest_file}")
    print(f"- {SUMMARY_FILE}")

    print("\nGenomic context summary:")
    print(context_summary.to_string(index=False))

    print("\nFirst-pass functional class summary:")
    print(class_summary.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()