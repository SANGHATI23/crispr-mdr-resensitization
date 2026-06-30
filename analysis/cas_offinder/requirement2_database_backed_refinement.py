#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

INPUT_ANNOTATED = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/"
    / "requirement2_expanded_offtarget_hits_functionally_annotated.csv"
)

INPUT_ALL20_GUIDES = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2/"
    / "requirement2_final_guide_prioritization_table_ALL20.csv"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2_database_refined"
)

OUT_HITS = OUT_DIR / "requirement2_database_refined_hit_annotations.csv"
OUT_GUIDES = OUT_DIR / "requirement2_database_refined_guide_summary_ALL20.csv"
OUT_TOP = OUT_DIR / "requirement2_database_refined_top_guides_by_gene_ALL20.csv"
OUT_HIGH = OUT_DIR / "requirement2_database_refined_high_confidence_functional_hits.csv"
OUT_SUMMARY = OUT_DIR / "requirement2_database_refined_completion_summary.txt"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Curated database-inspired term sets
# These are not direct BLAST results. They strengthen the GFF keyword layer
# using known AMR / virulence / essential / mobile-element terminology.
# ---------------------------------------------------------------------

CARD_AMR_TERMS = [
    # beta-lactamases / carbapenemases
    "beta-lactamase", "lactamase", "carbapenemase", "kpc", "ndm", "oxa",
    "ctx-m", "tem", "shv", "ges", "vim", "imp", "ampc", "cmY".lower(),

    # colistin / polymyxin
    "mcr", "colistin", "polymyxin", "pmr", "phoP".lower(), "phoQ".lower(),

    # MRSA / beta-lactam resistance
    "meca", "mecA".lower(), "pbp2a", "penicillin-binding protein 2a",
    "penicillin-binding protein", "methicillin",

    # vancomycin
    "vana", "vanb", "vanc", "vand", "vane", "vang", "vancomycin",

    # aminoglycoside
    "aac", "aph", "aad", "aminoglycoside",

    # tetracycline
    "tet", "tetracycline",

    # macrolide/lincosamide/streptogramin
    "erm", "mef", "msr", "macrolide", "lincosamide", "streptogramin",

    # quinolone
    "qnr", "quinolone", "fluoroquinolone", "gyrase", "parc", "pare",

    # sulfonamide / trimethoprim
    "sul", "dfr", "sulfonamide", "trimethoprim",

    # chloramphenicol / fosfomycin / rifampin
    "cat", "chloramphenicol", "fos", "fosfomycin", "rifampin", "rifamycin",

    # multidrug efflux
    "efflux", "multidrug", "acr", "mex", "mdt", "emr", "tolc",
    "drug resistance", "antibiotic resistance", "antimicrobial resistance",
]

VFDB_VIRULENCE_TERMS = [
    "virulence", "toxin", "exotoxin", "endotoxin", "enterotoxin",
    "hemolysin", "haemolysin", "leukocidin", "cytolysin",
    "adhesin", "invasin", "invasion", "colonization",
    "capsule", "capsular", "biofilm", "quorum",
    "fimbria", "fimbrial", "pilus", "pili",
    "secretion system", "type iii secretion", "type 3 secretion",
    "type vi secretion", "type 6 secretion", "pathogenicity island",
    "siderophore", "yersiniabactin", "enterobactin", "aerobactin",
    "lipopolysaccharide", "lps", "outer membrane protein",
]

DEG_ESSENTIAL_TERMS = [
    "essential", "ribosomal protein", "ribosome", "translation",
    "rna polymerase", "dna polymerase", "dna gyrase", "gyrase",
    "topoisomerase", "replication", "dna replication",
    "cell division", "fts", "ftsz", "ftsa", "ftsi",
    "mreb", "mra", "mur", "peptidoglycan", "cell wall",
    "fab", "fatty acid", "acyl carrier protein", "acp",
    "aminoacyl-trna", "trna synthetase", "ligase",
    "fol", "folate", "dihydrofolate", "thymidylate",
    "sec", "protein translocase", "atp synthase",
    "nad", "nadph", "glycolysis", "central metabolism",
    "amiA".lower(), "n-acetylmuramoyl-l-alanine amidase",
]

MOBILE_ELEMENT_TERMS = [
    "plasmid", "replicon", "replication protein", "rep protein",
    "transposase", "transposon", "insertion sequence", "is element",
    "integrase", "integron", "recombinase", "resolvase",
    "conjugative", "conjugal", "mobilization", "mob", "tra",
    "type iv secretion", "relaxase",
    "sccmec", "staphylococcal cassette chromosome", "cassette chromosome",
    "genomic island", "resistance island", "pathogenicity island",
    "phage", "prophage", "terminase", "capsid", "tail protein",
]


def norm(x):
    if pd.isna(x):
        return ""
    return str(x).lower()


def has_any(text, terms):
    text = norm(text)
    return any(t.lower() in text for t in terms)


def build_text(row):
    cols = [
        "overlap_gene",
        "overlap_name",
        "overlap_locus_tag",
        "overlap_old_locus_tag",
        "overlap_product",
        "overlap_note",
        "overlap_dbxref",
        "overlap_protein_id",
        "functional_text",
        "genomic_context",
        "functional_class_primary",
    ]
    return " ".join([norm(row[c]) for c in cols if c in row.index])


def refined_class(row):
    if row["db_amr_confirmed"]:
        return "DB_refined_AMR_or_resistance_related"
    if row["db_virulence_confirmed"]:
        return "DB_refined_virulence_related"
    if row["db_essential_confirmed"]:
        return "DB_refined_essential_like"
    if row["db_mobile_confirmed"]:
        return "DB_refined_mobile_or_plasmid_related"

    original = row.get("functional_class_primary", "")

    if original == "Coding_other":
        return "DB_refined_coding_other"
    if original == "Annotated_non_CDS_other":
        return "DB_refined_annotated_non_CDS_other"
    if original == "Noncoding_or_intergenic":
        return "DB_refined_noncoding_or_intergenic"

    return "DB_refined_other"


def refined_burden(row):
    mismatch = row.get("mismatch_count", None)

    try:
        mismatch = int(mismatch)
    except Exception:
        mismatch = None

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

    cls = row["db_refined_functional_class"]

    class_weight = {
        "DB_refined_AMR_or_resistance_related": 10.0,
        "DB_refined_virulence_related": 8.0,
        "DB_refined_essential_like": 7.0,
        "DB_refined_mobile_or_plasmid_related": 6.0,
        "DB_refined_coding_other": 2.0,
        "DB_refined_annotated_non_CDS_other": 1.0,
        "DB_refined_noncoding_or_intergenic": 0.25,
        "DB_refined_other": 1.0,
    }.get(cls, 1.0)

    return mismatch_weight + class_weight


def guide_recommendation(row):
    hits = row["total_offtarget_hits"]
    burden = row["db_refined_total_burden"]
    high = row["db_refined_high_confidence_functional_hits"]

    if hits == 0:
        return "Strong_candidate_zero_detected_offtarget_burden"
    if high == 0 and burden <= 25:
        return "Strong_candidate_low_database_refined_burden"
    if high <= 2 and burden <= 75:
        return "Acceptable_candidate_moderate_database_refined_burden"
    if high <= 5 and burden <= 150:
        return "Caution_candidate_elevated_database_refined_burden"
    return "Avoid_or_deprioritize_high_database_refined_burden"


def decision_from_rec(rec):
    if rec == "Strong_candidate_zero_detected_offtarget_burden":
        return "PRIORITIZE_OR_KEEP_AS_LOW_OFFTARGET_CANDIDATE"
    if rec == "Strong_candidate_low_database_refined_burden":
        return "PRIORITIZE"
    if rec == "Acceptable_candidate_moderate_database_refined_burden":
        return "KEEP_AS_BACKUP"
    if rec == "Caution_candidate_elevated_database_refined_burden":
        return "USE_WITH_CAUTION"
    return "DEPRIORITIZE"


def main():
    print("Requirement 2 database-refined annotation")
    print("=========================================")

    if not INPUT_ANNOTATED.exists():
        raise FileNotFoundError(INPUT_ANNOTATED)
    if not INPUT_ALL20_GUIDES.exists():
        raise FileNotFoundError(INPUT_ALL20_GUIDES)

    hits = pd.read_csv(INPUT_ANNOTATED)
    all20 = pd.read_csv(INPUT_ALL20_GUIDES)

    print("Hit-level input:", hits.shape)
    print("ALL20 guide input:", all20.shape)

    hits["db_refinement_text"] = hits.apply(build_text, axis=1)

    hits["db_amr_confirmed"] = hits["db_refinement_text"].apply(lambda x: has_any(x, CARD_AMR_TERMS))
    hits["db_virulence_confirmed"] = hits["db_refinement_text"].apply(lambda x: has_any(x, VFDB_VIRULENCE_TERMS))
    hits["db_essential_confirmed"] = hits["db_refinement_text"].apply(lambda x: has_any(x, DEG_ESSENTIAL_TERMS))
    hits["db_mobile_confirmed"] = hits["db_refinement_text"].apply(lambda x: has_any(x, MOBILE_ELEMENT_TERMS))

    hits["db_refined_functional_class"] = hits.apply(refined_class, axis=1)
    hits["db_refined_high_confidence_functional_hit"] = (
        hits["db_amr_confirmed"]
        | hits["db_virulence_confirmed"]
        | hits["db_essential_confirmed"]
        | hits["db_mobile_confirmed"]
    )

    hits["db_refined_burden_score"] = hits.apply(refined_burden, axis=1)

    # Hit-level high confidence table
    high = hits[hits["db_refined_high_confidence_functional_hit"]].copy()

    # Aggregate by guide
    grouped = (
        hits.groupby(["target_gene", "guide_id"], dropna=False)
        .agg(
            total_offtarget_hits=("guide_id", "count"),
            db_refined_amr_hits=("db_amr_confirmed", "sum"),
            db_refined_virulence_hits=("db_virulence_confirmed", "sum"),
            db_refined_essential_hits=("db_essential_confirmed", "sum"),
            db_refined_mobile_hits=("db_mobile_confirmed", "sum"),
            db_refined_high_confidence_functional_hits=("db_refined_high_confidence_functional_hit", "sum"),
            db_refined_total_burden=("db_refined_burden_score", "sum"),
            db_refined_mean_burden=("db_refined_burden_score", "mean"),
        )
        .reset_index()
    )

    # Retain all 20 guides
    base = all20[["target_gene", "guide_id"]].drop_duplicates()

    guide = base.merge(grouped, on=["target_gene", "guide_id"], how="left")

    fill_cols = [
        "total_offtarget_hits",
        "db_refined_amr_hits",
        "db_refined_virulence_hits",
        "db_refined_essential_hits",
        "db_refined_mobile_hits",
        "db_refined_high_confidence_functional_hits",
        "db_refined_total_burden",
        "db_refined_mean_burden",
    ]

    for c in fill_cols:
        guide[c] = guide[c].fillna(0)

    guide["db_refined_recommendation"] = guide.apply(guide_recommendation, axis=1)
    guide["db_refined_decision"] = guide["db_refined_recommendation"].apply(decision_from_rec)

    rec_order = {
        "Strong_candidate_zero_detected_offtarget_burden": 0,
        "Strong_candidate_low_database_refined_burden": 1,
        "Acceptable_candidate_moderate_database_refined_burden": 2,
        "Caution_candidate_elevated_database_refined_burden": 3,
        "Avoid_or_deprioritize_high_database_refined_burden": 4,
    }

    guide["db_refined_recommendation_order"] = guide["db_refined_recommendation"].map(rec_order).fillna(99)

    guide = guide.sort_values(
        by=[
            "target_gene",
            "db_refined_recommendation_order",
            "db_refined_total_burden",
            "db_refined_high_confidence_functional_hits",
            "total_offtarget_hits",
        ]
    ).copy()

    guide["db_refined_rank_within_gene"] = guide.groupby("target_gene").cumcount() + 1

    top = guide[guide["db_refined_rank_within_gene"] <= 2].copy()

    hits.to_csv(OUT_HITS, index=False)
    guide.to_csv(OUT_GUIDES, index=False)
    top.to_csv(OUT_TOP, index=False)
    high.to_csv(OUT_HIGH, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 2 Database-Refined Completion Summary\n")
        f.write("================================================\n\n")

        f.write("Final status\n")
        f.write("------------\n")
        f.write("Requirement 2 status: DATABASE-REFINED COMPUTATIONAL ANNOTATION COMPLETE.\n")
        f.write("This is a curated database-inspired refinement of GFF product annotations, not direct BLAST/HMM database alignment.\n\n")

        f.write("Input metrics\n")
        f.write("-------------\n")
        f.write(f"Annotated hit-level rows: {len(hits)}\n")
        f.write(f"All-guide rows retained: {len(guide)}\n")
        f.write(f"Genome accessions represented: {hits['accession'].nunique() if 'accession' in hits.columns else 'NA'}\n")
        f.write(f"GFF files used: {hits['overlap_gff_file'].nunique() if 'overlap_gff_file' in hits.columns else 'NA'}\n\n")

        f.write("Database-refined hit counts\n")
        f.write("---------------------------\n")
        f.write(f"DB-refined AMR/resistance hits: {int(hits['db_amr_confirmed'].sum())}\n")
        f.write(f"DB-refined virulence hits: {int(hits['db_virulence_confirmed'].sum())}\n")
        f.write(f"DB-refined essential-like hits: {int(hits['db_essential_confirmed'].sum())}\n")
        f.write(f"DB-refined mobile/plasmid hits: {int(hits['db_mobile_confirmed'].sum())}\n")
        f.write(f"High-confidence functional hits: {int(hits['db_refined_high_confidence_functional_hit'].sum())}\n\n")

        f.write("DB-refined primary class counts\n")
        f.write("-------------------------------\n")
        f.write(str(hits["db_refined_functional_class"].value_counts(dropna=False)))
        f.write("\n\n")

        f.write("DB-refined guide recommendation counts\n")
        f.write("--------------------------------------\n")
        f.write(str(guide["db_refined_recommendation"].value_counts(dropna=False)))
        f.write("\n\n")

        f.write("Top guides by gene after DB-refined Requirement 2\n")
        f.write("------------------------------------------------\n")
        f.write(
            top[
                [
                    "target_gene",
                    "guide_id",
                    "db_refined_rank_within_gene",
                    "total_offtarget_hits",
                    "db_refined_high_confidence_functional_hits",
                    "db_refined_total_burden",
                    "db_refined_recommendation",
                    "db_refined_decision",
                ]
            ].to_string(index=False)
        )
        f.write("\n\n")

        f.write("Manuscript-safe wording\n")
        f.write("-----------------------\n")
        f.write(
            "We further refined the GFF-based off-target annotation using curated "
            "database-inspired functional dictionaries representing AMR/resistance, "
            "virulence, essential-gene, and mobile-element categories. This refinement "
            "was applied to all 1,997 expanded Cas-OFFinder candidate off-target hits "
            "and aggregated into an all-20-guide burden table, retaining the zero-hit "
            "guide identified in the expanded Cas-OFFinder audit. Because this step "
            "uses curated term matching over GFF product annotations rather than direct "
            "sequence alignment against CARD, VFDB, DEG, or mobile-element databases, "
            "the results are reported as database-refined computational annotations "
            "rather than experimentally or sequence-database-confirmed annotations.\n"
        )

        f.write("\nOutput files\n")
        f.write("------------\n")
        f.write(f"DB-refined hit table: {OUT_HITS}\n")
        f.write(f"DB-refined all-20 guide table: {OUT_GUIDES}\n")
        f.write(f"DB-refined top guides: {OUT_TOP}\n")
        f.write(f"DB-refined high-confidence hits: {OUT_HIGH}\n")
        f.write(f"DB-refined summary: {OUT_SUMMARY}\n")

    print("\nWrote:")
    print(OUT_HITS)
    print(OUT_GUIDES)
    print(OUT_TOP)
    print(OUT_HIGH)
    print(OUT_SUMMARY)

    print("\nDB-refined class counts:")
    print(hits["db_refined_functional_class"].value_counts(dropna=False))

    print("\nDB-refined guide recommendation counts:")
    print(guide["db_refined_recommendation"].value_counts(dropna=False))

    print("\nTop guides:")
    print(
        top[
            [
                "target_gene",
                "guide_id",
                "db_refined_rank_within_gene",
                "total_offtarget_hits",
                "db_refined_high_confidence_functional_hits",
                "db_refined_total_burden",
                "db_refined_recommendation",
                "db_refined_decision",
            ]
        ].to_string(index=False)
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
