#!/usr/bin/env python3

"""
Requirement 2:
Annotate expanded Cas-OFFinder off-target hits using matched GFF/GFF3 files
and keyword-based functional classes.

Input:
results_cas_offinder/expanded_panel/final_requirement1/
    requirement1_expanded_casoffinder_genomewide_offtarget_hit_table.csv

GFF/GFF3 folder:
data_cas_offinder/expanded_panel/annotations/

Outputs:
results_cas_offinder/expanded_panel/functional_annotation/
    requirement2_expanded_offtarget_hits_functionally_annotated.csv
    requirement2_high_interest_functional_offtargets.csv
    requirement2_guide_level_functional_offtarget_burden.csv
    requirement2_functional_annotation_summary.txt
"""

from pathlib import Path
import pandas as pd
import re


PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

INPUT_HITS = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/final_requirement1/"
    / "requirement1_expanded_casoffinder_genomewide_offtarget_hit_table.csv"
)

GFF_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/annotations"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/functional_annotation"

OUT_ANNOTATED = OUT_DIR / "requirement2_expanded_offtarget_hits_functionally_annotated.csv"
OUT_HIGH_INTEREST = OUT_DIR / "requirement2_high_interest_functional_offtargets.csv"
OUT_GUIDE_SUMMARY = OUT_DIR / "requirement2_guide_level_functional_offtarget_burden.csv"
OUT_SUMMARY = OUT_DIR / "requirement2_functional_annotation_summary.txt"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------
# Keyword dictionaries
# -----------------------------

AMR_KEYWORDS = [
    "beta-lactamase",
    "lactamase",
    "carbapenemase",
    "kpc",
    "ndm",
    "mcr",
    "meca",
    "pbp2a",
    "penicillin-binding protein 2a",
    "vancomycin",
    "vana",
    "vanb",
    "resistance",
    "resistant",
    "antibiotic",
    "antimicrobial",
    "efflux",
    "multidrug",
    "drug resistance",
    "aminoglycoside",
    "tetracycline",
    "macrolide",
    "quinolone",
    "fluoroquinolone",
    "sulfonamide",
    "chloramphenicol",
    "fosfomycin",
    "colistin",
    "polymyxin",
]

VIRULENCE_KEYWORDS = [
    "virulence",
    "toxin",
    "hemolysin",
    "haemolysin",
    "adhesin",
    "invasion",
    "invasin",
    "capsule",
    "capsular",
    "biofilm",
    "fimbrial",
    "fimbria",
    "pilus",
    "pili",
    "secretion system",
    "type iii secretion",
    "type 3 secretion",
    "type vi secretion",
    "type 6 secretion",
    "pathogenicity",
    "enterotoxin",
    "exotoxin",
    "leukocidin",
    "siderophore",
]

ESSENTIAL_KEYWORDS = [
    "dna polymerase",
    "rna polymerase",
    "ribosomal protein",
    "ribosome",
    "gyrase",
    "topoisomerase",
    "cell division",
    "fts",
    "mreb",
    "peptidoglycan",
    "mur",
    "ddl",
    "fab",
    "acp",
    "dna replication",
    "replication initiation",
    "protein synthesis",
    "translation",
    "transcription",
    "essential",
    "trna synthetase",
    "aminoacyl-trna",
]

MOBILE_KEYWORDS = [
    "transposase",
    "integrase",
    "recombinase",
    "insertion sequence",
    "is element",
    "mobile element",
    "plasmid",
    "conjugative",
    "mobilization",
    "mobility",
    "resistance island",
    "genomic island",
    "integron",
    "sccmec",
    "cassette chromosome",
    "phage",
    "prophage",
    "replicon",
]


# -----------------------------
# Helper functions
# -----------------------------

def normalize_text(x):
    if pd.isna(x):
        return ""
    return str(x).lower()


def contains_keyword(text, keywords):
    text = normalize_text(text)
    return any(k.lower() in text for k in keywords)


def parse_gff_attributes(attr_text):
    attrs = {}

    if pd.isna(attr_text):
        return attrs

    attr_text = str(attr_text)

    for item in attr_text.split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')

    return attrs


def load_gff3_files(gff_dir):
    print(f"Loading GFF/GFF3 files from: {gff_dir}")

    gff_files = sorted(list(gff_dir.glob("*.gff")) + list(gff_dir.glob("*.gff3")))

    if not gff_files:
        raise FileNotFoundError(f"No .gff or .gff3 files found in: {gff_dir}")

    all_records = []

    for gff_file in gff_files:
        assembly_accession = gff_file.name.replace("_genomic.gff", "").replace("_genomic.gff3", "")

        with open(gff_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")

                if len(parts) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    continue

                attrs = parse_gff_attributes(attributes)

                gene = attrs.get("gene", "")
                name = attrs.get("Name", "")
                locus_tag = attrs.get("locus_tag", "")
                product = attrs.get("product", "")
                note = attrs.get("Note", "")
                dbxref = attrs.get("Dbxref", "")
                protein_id = attrs.get("protein_id", "")
                old_locus_tag = attrs.get("old_locus_tag", "")

                searchable_text = " ".join(
                    [
                        gene,
                        name,
                        locus_tag,
                        old_locus_tag,
                        product,
                        note,
                        dbxref,
                        protein_id,
                        feature_type,
                    ]
                ).lower()

                all_records.append(
                    {
                        "gff_file": gff_file.name,
                        "gff_assembly_accession": assembly_accession,
                        "seqid": seqid,
                        "feature_type": feature_type,
                        "feature_start": start,
                        "feature_end": end,
                        "feature_strand": strand,
                        "gene": gene,
                        "name": name,
                        "locus_tag": locus_tag,
                        "old_locus_tag": old_locus_tag,
                        "product": product,
                        "note": note,
                        "dbxref": dbxref,
                        "protein_id": protein_id,
                        "searchable_text": searchable_text,
                    }
                )

    gff_df = pd.DataFrame(all_records)

    print(f"Loaded GFF feature rows: {len(gff_df)}")
    print(f"Loaded GFF files: {len(gff_files)}")

    return gff_df


def infer_column(df, candidates, required=True):
    for col in candidates:
        if col in df.columns:
            return col

    if required:
        raise ValueError(
            "Could not infer required column. Tried: "
            + ", ".join(candidates)
            + "\nAvailable columns are:\n"
            + "\n".join(df.columns)
        )

    return None


def infer_hit_seqid(row):
    possible_cols = [
        "seqid",
        "contig",
        "contig_id",
        "contig_description",
        "chromosome",
        "genome_contig",
        "reference_sequence",
        "target_contig",
        "subject_id",
        "reference_id",
    ]

    for col in possible_cols:
        if col in row.index and pd.notna(row[col]):
            value = str(row[col]).strip()
            if value:
                return value.split()[0]

    return ""


def infer_hit_accession(row):
    possible_cols = [
        "accession",
        "genome_accession",
        "assembly_accession",
        "genome_id",
        "genome_label",
        "assembly",
    ]

    for col in possible_cols:
        if col in row.index and pd.notna(row[col]):
            value = str(row[col]).strip()
            if value:
                # Pull GCF/GCA style accession if present
                match = re.search(r"(GC[AF]_\d+\.\d+)", value)
                if match:
                    return match.group(1)
                return value

    # Fallback from contig description
    for col in ["contig_description", "contig", "seqid"]:
        if col in row.index and pd.notna(row[col]):
            value = str(row[col])
            match = re.search(r"(GC[AF]_\d+\.\d+)", value)
            if match:
                return match.group(1)

    return ""


def infer_hit_position(row):
    possible_cols = [
        "position",
        "position_1based",
        "genomic_position",
        "hit_position",
        "start",
        "offtarget_position",
        "casoffinder_position",
    ]

    for col in possible_cols:
        if col in row.index and pd.notna(row[col]):
            try:
                return int(row[col])
            except Exception:
                continue

    raise ValueError(
        "Could not infer hit position column.\nAvailable columns:\n"
        + "\n".join(row.index)
    )


def infer_mismatch(row):
    possible_cols = [
        "mismatches",
        "mismatch_count",
        "n_mismatches",
        "mm",
        "mismatch",
    ]

    for col in possible_cols:
        if col in row.index and pd.notna(row[col]):
            try:
                return int(row[col])
            except Exception:
                continue

    return None


def choose_best_overlap(overlap_df):
    """
    If multiple features overlap, prefer CDS, then gene, then first feature.
    """
    if overlap_df.empty:
        return None

    cds = overlap_df[overlap_df["feature_type"].str.lower() == "cds"]
    if not cds.empty:
        return cds.iloc[0]

    gene = overlap_df[overlap_df["feature_type"].str.lower() == "gene"]
    if not gene.empty:
        return gene.iloc[0]

    return overlap_df.iloc[0]


def annotate_hits_with_gff(hits_df, gff_df):
    print("Indexing GFF features...")

    gff_by_seqid = {
        seqid: sub.copy()
        for seqid, sub in gff_df.groupby("seqid")
    }

    gff_by_accession = {
        acc: sub.copy()
        for acc, sub in gff_df.groupby("gff_assembly_accession")
    }

    annotated_rows = []

    print("Annotating off-target hits...")

    for i, row in hits_df.iterrows():
        if (i + 1) % 250 == 0:
            print(f"Annotated {i + 1} / {len(hits_df)} hits...")

        hit_seqid = infer_hit_seqid(row)
        hit_accession = infer_hit_accession(row)
        hit_pos = infer_hit_position(row)

        candidates = pd.DataFrame()

        # Best matching path: exact contig/seqid
        if hit_seqid and hit_seqid in gff_by_seqid:
            candidates = gff_by_seqid[hit_seqid]

        # Fallback: assembly accession
        if candidates.empty and hit_accession in gff_by_accession:
            candidates = gff_by_accession[hit_accession]

        # Final fallback: file name contains accession
        if candidates.empty and hit_accession:
            candidates = gff_df[
                gff_df["gff_file"].str.contains(hit_accession, case=False, na=False)
            ]

        if not candidates.empty:
            overlap = candidates[
                (candidates["feature_start"] <= hit_pos)
                & (candidates["feature_end"] >= hit_pos)
            ]
        else:
            overlap = pd.DataFrame()

        out = row.to_dict()

        out["annotation_hit_seqid_used"] = hit_seqid
        out["annotation_hit_accession_used"] = hit_accession
        out["annotation_hit_position_used"] = hit_pos

        if overlap.empty:
            out["genomic_context"] = "Noncoding_or_intergenic"
            out["overlap_feature_type"] = ""
            out["overlap_gene"] = ""
            out["overlap_name"] = ""
            out["overlap_locus_tag"] = ""
            out["overlap_old_locus_tag"] = ""
            out["overlap_product"] = ""
            out["overlap_note"] = ""
            out["overlap_dbxref"] = ""
            out["overlap_protein_id"] = ""
            out["overlap_feature_start"] = ""
            out["overlap_feature_end"] = ""
            out["overlap_feature_strand"] = ""
            out["overlap_gff_file"] = ""
            out["functional_text"] = ""

        else:
            best = choose_best_overlap(overlap)

            feature_type = str(best.get("feature_type", ""))

            if feature_type.lower() == "cds":
                genomic_context = "Coding_CDS"
            elif feature_type.lower() == "gene":
                genomic_context = "Annotated_gene_feature"
            else:
                genomic_context = "Annotated_non_CDS_feature"

            functional_text = " ".join(
                [
                    str(best.get("gene", "")),
                    str(best.get("name", "")),
                    str(best.get("locus_tag", "")),
                    str(best.get("old_locus_tag", "")),
                    str(best.get("product", "")),
                    str(best.get("note", "")),
                    str(best.get("dbxref", "")),
                    str(best.get("protein_id", "")),
                    str(best.get("feature_type", "")),
                ]
            )

            out["genomic_context"] = genomic_context
            out["overlap_feature_type"] = best.get("feature_type", "")
            out["overlap_gene"] = best.get("gene", "")
            out["overlap_name"] = best.get("name", "")
            out["overlap_locus_tag"] = best.get("locus_tag", "")
            out["overlap_old_locus_tag"] = best.get("old_locus_tag", "")
            out["overlap_product"] = best.get("product", "")
            out["overlap_note"] = best.get("note", "")
            out["overlap_dbxref"] = best.get("dbxref", "")
            out["overlap_protein_id"] = best.get("protein_id", "")
            out["overlap_feature_start"] = best.get("feature_start", "")
            out["overlap_feature_end"] = best.get("feature_end", "")
            out["overlap_feature_strand"] = best.get("feature_strand", "")
            out["overlap_gff_file"] = best.get("gff_file", "")
            out["functional_text"] = functional_text

        text = normalize_text(out.get("functional_text", ""))

        out["is_amr_keyword"] = contains_keyword(text, AMR_KEYWORDS)
        out["is_virulence_keyword"] = contains_keyword(text, VIRULENCE_KEYWORDS)
        out["is_essential_keyword"] = contains_keyword(text, ESSENTIAL_KEYWORDS)
        out["is_mobile_keyword"] = contains_keyword(text, MOBILE_KEYWORDS)

        if out["is_amr_keyword"]:
            primary_class = "AMR_or_resistance_related"
        elif out["is_virulence_keyword"]:
            primary_class = "Virulence_related"
        elif out["is_essential_keyword"]:
            primary_class = "Essential_like"
        elif out["is_mobile_keyword"]:
            primary_class = "Mobile_or_plasmid_related"
        elif out["genomic_context"] == "Noncoding_or_intergenic":
            primary_class = "Noncoding_or_intergenic"
        elif out["genomic_context"] == "Coding_CDS":
            primary_class = "Coding_other"
        elif out["genomic_context"] == "Annotated_gene_feature":
            primary_class = "Annotated_gene_other"
        else:
            primary_class = "Annotated_non_CDS_other"

        out["functional_class_primary"] = primary_class

        mismatch = infer_mismatch(row)

        out["is_low_mismatch_hit"] = False if mismatch is None else mismatch <= 2

        out["is_high_interest_functional_hit"] = (
            primary_class
            in [
                "AMR_or_resistance_related",
                "Virulence_related",
                "Essential_like",
                "Mobile_or_plasmid_related",
            ]
            or out["is_low_mismatch_hit"]
        )

        annotated_rows.append(out)

    return pd.DataFrame(annotated_rows)


def compute_functional_burden(row):
    mismatch = infer_mismatch(row)

    if mismatch is None:
        mismatch_weight = 1.0
    elif mismatch <= 1:
        mismatch_weight = 5.0
    elif mismatch == 2:
        mismatch_weight = 3.0
    elif mismatch == 3:
        mismatch_weight = 2.0
    else:
        mismatch_weight = 1.0

    functional_class = row.get("functional_class_primary", "")

    class_weight = {
        "AMR_or_resistance_related": 8.0,
        "Virulence_related": 7.0,
        "Essential_like": 6.0,
        "Mobile_or_plasmid_related": 5.0,
        "Coding_other": 2.0,
        "Annotated_gene_other": 1.5,
        "Annotated_non_CDS_other": 1.0,
        "Noncoding_or_intergenic": 0.25,
    }.get(functional_class, 1.0)

    return mismatch_weight + class_weight


def infer_guide_column(df):
    candidates = [
        "guide_id",
        "selected_guide_id",
        "guide_name",
        "query_id",
        "guide",
    ]

    return infer_column(df, candidates, required=True)


def infer_gene_column(df):
    candidates = [
        "target_gene",
        "gene",
        "amr_gene",
        "target",
    ]

    return infer_column(df, candidates, required=False)


def make_guide_summary(annotated):
    guide_col = infer_guide_column(annotated)
    gene_col = infer_gene_column(annotated)

    group_cols = [guide_col]
    if gene_col:
        group_cols = [gene_col, guide_col]

    guide_summary = (
        annotated.groupby(group_cols, dropna=False)
        .agg(
            total_offtarget_hits=("functional_class_primary", "count"),
            coding_cds_hits=("genomic_context", lambda x: (x == "Coding_CDS").sum()),
            noncoding_or_intergenic_hits=(
                "genomic_context",
                lambda x: (x == "Noncoding_or_intergenic").sum(),
            ),
            amr_hits=("is_amr_keyword", "sum"),
            virulence_hits=("is_virulence_keyword", "sum"),
            essential_like_hits=("is_essential_keyword", "sum"),
            mobile_or_plasmid_hits=("is_mobile_keyword", "sum"),
            low_mismatch_hits=("is_low_mismatch_hit", "sum"),
            high_interest_hits=("is_high_interest_functional_hit", "sum"),
            total_functional_burden=("functional_burden_score", "sum"),
            mean_functional_burden=("functional_burden_score", "mean"),
        )
        .reset_index()
    )

    def recommendation(row):
        burden = row["total_functional_burden"]
        high = row["high_interest_hits"]

        if high == 0 and burden <= 25:
            return "Strong_candidate_low_functional_offtarget_burden"
        elif high <= 2 and burden <= 75:
            return "Acceptable_candidate_moderate_functional_offtarget_burden"
        elif high <= 5 and burden <= 150:
            return "Caution_candidate_elevated_functional_offtarget_burden"
        else:
            return "Avoid_or_deprioritize_high_functional_offtarget_burden"

    guide_summary["functional_annotation_recommendation"] = guide_summary.apply(
        recommendation, axis=1
    )

    sort_cols = []
    if gene_col:
        sort_cols.append(gene_col)

    sort_cols.extend(
        [
            "functional_annotation_recommendation",
            "total_functional_burden",
            "high_interest_hits",
            "total_offtarget_hits",
        ]
    )

    guide_summary = guide_summary.sort_values(
        by=sort_cols,
        ascending=True,
    )

    return guide_summary


def write_summary(hits, gff, annotated, high_interest, guide_summary):
    total_hits = len(annotated)
    labeled_hits = annotated["functional_class_primary"].notna().sum()
    blank_labels = annotated["functional_class_primary"].isna().sum()

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 2 Functional Annotation Summary\n")
        f.write("===========================================\n\n")

        f.write("Purpose\n")
        f.write("-------\n")
        f.write(
            "This step annotates the expanded Cas-OFFinder off-target hit table "
            "using matched GFF/GFF3 genome annotation files and keyword-derived "
            "functional classes.\n\n"
        )

        f.write("Input and annotation files\n")
        f.write("--------------------------\n")
        f.write(f"Input expanded off-target hit table: {INPUT_HITS}\n")
        f.write(f"GFF/GFF3 annotation folder: {GFF_DIR}\n")
        f.write(f"Input expanded off-target hits: {len(hits)}\n")
        f.write(f"GFF feature rows loaded: {len(gff)}\n")
        f.write(f"GFF files represented: {gff['gff_file'].nunique()}\n\n")

        f.write("Output QC\n")
        f.write("---------\n")
        f.write(f"Annotated output rows: {len(annotated)}\n")
        f.write(f"Rows with functional label: {labeled_hits}\n")
        f.write(f"Rows missing functional label: {blank_labels}\n")
        f.write(f"High-interest functional hits: {len(high_interest)}\n")
        f.write(f"Guide-level rows: {len(guide_summary)}\n\n")

        f.write("Genomic context counts\n")
        f.write("----------------------\n")
        f.write(str(annotated["genomic_context"].value_counts(dropna=False)))
        f.write("\n\n")

        f.write("Primary functional class counts\n")
        f.write("-------------------------------\n")
        f.write(str(annotated["functional_class_primary"].value_counts(dropna=False)))
        f.write("\n\n")

        f.write("Keyword class counts\n")
        f.write("--------------------\n")
        f.write(f"AMR/resistance keyword hits: {int(annotated['is_amr_keyword'].sum())}\n")
        f.write(f"Virulence keyword hits: {int(annotated['is_virulence_keyword'].sum())}\n")
        f.write(f"Essential-like keyword hits: {int(annotated['is_essential_keyword'].sum())}\n")
        f.write(f"Mobile/plasmid keyword hits: {int(annotated['is_mobile_keyword'].sum())}\n")
        f.write(f"Low-mismatch hits <=2 mismatches: {int(annotated['is_low_mismatch_hit'].sum())}\n\n")

        f.write("Guide-level recommendation counts\n")
        f.write("---------------------------------\n")
        f.write(
            str(
                guide_summary[
                    "functional_annotation_recommendation"
                ].value_counts(dropna=False)
            )
        )
        f.write("\n\n")

        f.write("Output files\n")
        f.write("------------\n")
        f.write(f"Annotated hit-level table: {OUT_ANNOTATED}\n")
        f.write(f"High-interest functional hits: {OUT_HIGH_INTEREST}\n")
        f.write(f"Guide-level functional burden table: {OUT_GUIDE_SUMMARY}\n")
        f.write(f"Functional annotation summary: {OUT_SUMMARY}\n\n")

        if blank_labels == 0 and total_hits == 1997:
            f.write(
                "Requirement 2 first-pass status: PASS for GFF3 + keyword-based "
                "annotation of all 1,997 expanded Cas-OFFinder hits.\n"
            )
        elif blank_labels == 0:
            f.write(
                "Requirement 2 first-pass status: PASS for all rows present in the input file. "
                "Note: input row count was not exactly 1,997, so confirm the input file.\n"
            )
        else:
            f.write(
                "Requirement 2 first-pass status: CHECK REQUIRED because some rows are missing labels.\n"
            )

        f.write("\nManuscript-safe wording\n")
        f.write("-----------------------\n")
        f.write(
            "We extended the expanded Cas-OFFinder hit table with matched GFF3-based "
            "functional annotation. Each candidate off-target hit was assigned a genomic "
            "context and first-pass functional class using overlapping genome features and "
            "keyword-derived biological categories, including AMR/resistance-related, "
            "virulence-related, essential-like, mobile/plasmid-related, coding-other, and "
            "noncoding/intergenic classes. These annotations were aggregated into a "
            "guide-level functional off-target burden table for downstream FOTR-CRISPR "
            "prioritization. This annotation remains a first-pass computational "
            "classification and should be strengthened with database-backed CARD, VFDB, "
            "DEG, and mobile-element mapping in the next stage.\n"
        )


def main():
    print("Requirement 2: Expanded off-target functional annotation")
    print("=======================================================")

    if not INPUT_HITS.exists():
        raise FileNotFoundError(f"Input hit table not found: {INPUT_HITS}")

    if not GFF_DIR.exists():
        raise FileNotFoundError(f"GFF annotation folder not found: {GFF_DIR}")

    print(f"Input hit table: {INPUT_HITS}")
    hits = pd.read_csv(INPUT_HITS)
    print(f"Input hits shape: {hits.shape}")

    print("\nHit table columns:")
    for col in hits.columns:
        print(f"- {col}")

    gff = load_gff3_files(GFF_DIR)

    annotated = annotate_hits_with_gff(hits, gff)

    print("\nCalculating functional burden scores...")
    annotated["functional_burden_score"] = annotated.apply(
        compute_functional_burden, axis=1
    )

    high_interest = annotated[
        annotated["is_high_interest_functional_hit"] == True
    ].copy()

    guide_summary = make_guide_summary(annotated)

    annotated.to_csv(OUT_ANNOTATED, index=False)
    high_interest.to_csv(OUT_HIGH_INTEREST, index=False)
    guide_summary.to_csv(OUT_GUIDE_SUMMARY, index=False)
    write_summary(hits, gff, annotated, high_interest, guide_summary)

    print("\nWrote output files:")
    print(f"- {OUT_ANNOTATED}")
    print(f"- {OUT_HIGH_INTEREST}")
    print(f"- {OUT_GUIDE_SUMMARY}")
    print(f"- {OUT_SUMMARY}")

    print("\nGenomic context counts:")
    print(annotated["genomic_context"].value_counts(dropna=False))

    print("\nPrimary functional class counts:")
    print(annotated["functional_class_primary"].value_counts(dropna=False))

    print("\nGuide recommendation counts:")
    print(guide_summary["functional_annotation_recommendation"].value_counts(dropna=False))

    print("\nDone.")


if __name__ == "__main__":
    main()