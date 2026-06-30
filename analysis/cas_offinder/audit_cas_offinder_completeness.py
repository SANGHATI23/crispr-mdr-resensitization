#!/usr/bin/env python3
"""
Audit completeness of the first-pass Cas-OFFinder analysis.

This script checks:

1. Whether all 20 selected guides are represented.
2. Which guides have zero detected off-target hits.
3. Whether all parsed hits map to a valid guide.
4. Whether mismatch counts are within the configured threshold.
5. Whether all expected genomes are represented.
6. Whether duplicate parsed rows exist.
7. Hit counts by guide, gene, genome, mismatch count, and strand.
8. Whether guide-level summary files accidentally omit zero-hit guides.

Outputs:
    results_cas_offinder/qc/cas_offinder_guide_completeness_audit.csv
    results_cas_offinder/qc/cas_offinder_zero_hit_guides.csv
    results_cas_offinder/qc/cas_offinder_duplicate_hits.csv
    results_cas_offinder/qc/cas_offinder_unmapped_hits.csv
    results_cas_offinder/qc/cas_offinder_all_hits_audited.csv
    results_cas_offinder/qc/cas_offinder_qc_summary.txt
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

SELECTED_GUIDES_FILE = (
    PROJECT_ROOT
    / "results_final_selection"
    / "top5_per_gene_risk_aware_guides.csv"
)

PARSED_HITS_FILE = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "parsed"
    / "cas_offinder_all_genome_offtarget_hits_parsed.csv"
)

EXISTING_GUIDE_SUMMARY_FILE = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "guide_level_fotr"
    / "guide_level_functional_offtarget_risk_summary.csv"
)

OUTPUT_DIR = PROJECT_ROOT / "results_cas_offinder" / "qc"

EXPECTED_GENES = {
    "blaKPC",
    "blaNDM1",
    "mcr1",
    "mecA",
}

EXPECTED_GENOMES = {
    "Klebsiella_pneumoniae_reference",
    "Ecoli_reference",
    "Staphylococcus_aureus_MRSA_reference",
}

MAX_ALLOWED_MISMATCHES = 4
EXPECTED_GUIDE_COUNT = 20


def normalize_column_name(column: str) -> str:
    """Convert a column name to a normalized comparison form."""
    return (
        str(column)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
        .replace("/", "_")
    )


def find_column(
    dataframe: pd.DataFrame,
    candidates: list[str],
    required: bool = True,
) -> str | None:
    """
    Find a dataframe column using several possible column-name candidates.
    """
    normalized_lookup = {
        normalize_column_name(column): column
        for column in dataframe.columns
    }

    for candidate in candidates:
        normalized_candidate = normalize_column_name(candidate)

        if normalized_candidate in normalized_lookup:
            return normalized_lookup[normalized_candidate]

    if required:
        raise KeyError(
            "Could not find any of these columns:\n"
            f"{candidates}\n\n"
            f"Available columns:\n{list(dataframe.columns)}"
        )

    return None


def clean_sequence(series: pd.Series) -> pd.Series:
    """Normalize DNA sequence values."""
    return (
        series.astype(str)
        .str.upper()
        .str.strip()
        .str.replace(r"\s+", "", regex=True)
    )


def clean_text(series: pd.Series) -> pd.Series:
    """Normalize general text values."""
    return series.astype(str).str.strip()


def fail(message: str) -> None:
    """Exit with a readable error."""
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    print("=" * 72)
    print("Cas-OFFinder completeness audit")
    print("=" * 72)
    print(f"Project root: {PROJECT_ROOT}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not SELECTED_GUIDES_FILE.exists():
        fail(
            "Selected-guide file not found:\n"
            f"{SELECTED_GUIDES_FILE}"
        )

    if not PARSED_HITS_FILE.exists():
        fail(
            "Parsed Cas-OFFinder file not found:\n"
            f"{PARSED_HITS_FILE}"
        )

    guides = pd.read_csv(SELECTED_GUIDES_FILE)
    hits = pd.read_csv(PARSED_HITS_FILE)

    print(f"\nLoaded selected guides: {guides.shape}")
    print(f"Loaded parsed hits:     {hits.shape}")

    print("\nSelected-guide columns:")
    print(list(guides.columns))

    print("\nParsed-hit columns:")
    print(list(hits.columns))

    # ---------------------------------------------------------
    # Identify selected-guide columns
    # ---------------------------------------------------------
    guide_id_col = find_column(
        guides,
        [
            "final_guide_id",
            "guide_id",
            "riskaware_guide_id",
            "selected_guide_id",
            "candidate_id",
        ],
    )

    guide_gene_col = find_column(
        guides,
        [
            "mlcb_gene",
            "gene",
            "target_gene",
            "amr_gene",
            "gene_name",
        ],
    )

    guide_spacer_col = find_column(
        guides,
        [
            "mlcb_spacer",
            "spacer",
            "guide_sequence",
            "spacer_sequence",
            "guide",
            "sequence",
        ],
    )

    guide_rank_col = find_column(
        guides,
        [
            "final_gene_rank",
            "risk_aware_rank",
            "rank",
            "guide_rank",
            "within_gene_rank",
            "functional_risk_rank",
        ],
        required=False,
    )

    # ---------------------------------------------------------
    # Identify parsed-hit columns
    # ---------------------------------------------------------
    hit_guide_id_col = find_column(
        hits,
        [
            "final_guide_id",
            "guide_id",
            "riskaware_guide_id",
            "selected_guide_id",
            "candidate_id",
        ],
        required=False,
    )

    hit_gene_col = find_column(
        hits,
        [
            "mlcb_gene",
            "gene",
            "target_gene",
            "amr_gene",
            "gene_name",
        ],
        required=False,
    )

    hit_spacer_col = find_column(
        hits,
        [
            "query_spacer",
            "mlcb_spacer",
            "spacer",
            "guide_sequence",
            "spacer_sequence",
            "query_sequence",
            "query",
        ],
    )

    hit_genome_col = find_column(
        hits,
        [
            "genome_label",
            "genome",
            "reference_genome",
            "genome_name",
            "source_genome",
        ],
    )

    hit_contig_col = find_column(
        hits,
        [
            "contig_description",
            "contig",
            "chromosome",
            "seqid",
            "reference_sequence",
        ],
        required=False,
    )

    hit_position_col = find_column(
        hits,
        [
            "position",
            "position_0based",
            "position_1based",
            "genomic_position",
            "start",
        ],
    )

    hit_sequence_col = find_column(
        hits,
        [
            "off_target_sequence",
            "offtarget_sequence",
            "matched_sequence",
            "dna_sequence",
            "target_sequence",
        ],
    )

    hit_strand_col = find_column(
        hits,
        [
            "strand",
            "target_strand",
        ],
    )

    hit_mismatch_col = find_column(
        hits,
        [
            "mismatch_count",
            "mismatches",
            "number_of_mismatches",
            "n_mismatches",
        ],
    )

    # ---------------------------------------------------------
    # Normalize selected-guide table
    # ---------------------------------------------------------
    guides_clean = pd.DataFrame(
        {
            "guide_id": clean_text(guides[guide_id_col]),
            "gene": clean_text(guides[guide_gene_col]),
            "spacer": clean_sequence(guides[guide_spacer_col]),
        }
    )

    if guide_rank_col is not None:
        guides_clean["original_rank"] = pd.to_numeric(
            guides[guide_rank_col],
            errors="coerce",
        )
    else:
        guides_clean["original_rank"] = pd.NA

    duplicate_selected_guides = guides_clean[
        guides_clean.duplicated(
            subset=["gene", "spacer"],
            keep=False,
        )
    ].sort_values(
        ["gene", "spacer"]
    )

    if not duplicate_selected_guides.empty:
        print(
            "\nWARNING: Duplicate gene-spacer rows exist "
            "in the selected guides:"
        )
        print(
            duplicate_selected_guides.to_string(
                index=False
            )
        )

    guides_clean = guides_clean.drop_duplicates(
        subset=["gene", "spacer"],
        keep="first",
    ).copy()

    # ---------------------------------------------------------
    # Normalize parsed-hit table
    # ---------------------------------------------------------
    hits_clean = hits.copy()

    hits_clean["spacer_normalized"] = clean_sequence(
        hits_clean[hit_spacer_col]
    )

    # Cas-OFFinder query values may contain spacer + NNN.
    hits_clean["spacer_normalized"] = (
        hits_clean["spacer_normalized"]
        .str.replace(r"NNN$", "", regex=True)
        .str[:20]
    )

    hits_clean["genome_normalized"] = clean_text(
        hits_clean[hit_genome_col]
    )

    hits_clean["offtarget_sequence_normalized"] = clean_sequence(
        hits_clean[hit_sequence_col]
    )

    hits_clean["strand_normalized"] = clean_text(
        hits_clean[hit_strand_col]
    )

    hits_clean["position_normalized"] = pd.to_numeric(
        hits_clean[hit_position_col],
        errors="coerce",
    )

    hits_clean["mismatch_count_normalized"] = pd.to_numeric(
        hits_clean[hit_mismatch_col],
        errors="coerce",
    )

    if hit_guide_id_col is not None:
        hits_clean["existing_guide_id"] = clean_text(
            hits_clean[hit_guide_id_col]
        )
    else:
        hits_clean["existing_guide_id"] = pd.NA

    if hit_gene_col is not None:
        hits_clean["existing_gene"] = clean_text(
            hits_clean[hit_gene_col]
        )
    else:
        hits_clean["existing_gene"] = pd.NA

    # ---------------------------------------------------------
    # Authoritative guide mapping using selected spacer table
    # ---------------------------------------------------------
    spacer_map = guides_clean[
        [
            "guide_id",
            "gene",
            "spacer",
            "original_rank",
        ]
    ].rename(
        columns={
            "spacer": "spacer_normalized",
        }
    )

    audited_hits = hits_clean.merge(
        spacer_map,
        on="spacer_normalized",
        how="left",
        validate="many_to_one",
    )

    # ---------------------------------------------------------
    # Mapping checks
    # ---------------------------------------------------------
    unmapped_hits = audited_hits[
        audited_hits["guide_id"].isna()
    ].copy()

    conflicting_guide_ids = audited_hits[
        audited_hits["existing_guide_id"].notna()
        & audited_hits["guide_id"].notna()
        & (
            audited_hits["existing_guide_id"]
            != audited_hits["guide_id"]
        )
    ].copy()

    conflicting_genes = audited_hits[
        audited_hits["existing_gene"].notna()
        & audited_hits["gene"].notna()
        & (
            audited_hits["existing_gene"]
            != audited_hits["gene"]
        )
    ].copy()

    # ---------------------------------------------------------
    # Duplicate hit checks
    # ---------------------------------------------------------
    duplicate_key = [
        "guide_id",
        "genome_normalized",
        "position_normalized",
        "offtarget_sequence_normalized",
        "strand_normalized",
        "mismatch_count_normalized",
    ]

    if hit_contig_col is not None:
        audited_hits["contig_normalized"] = clean_text(
            audited_hits[hit_contig_col]
        )

        duplicate_key.insert(
            2,
            "contig_normalized",
        )

    duplicate_hits = audited_hits[
        audited_hits.duplicated(
            subset=duplicate_key,
            keep=False,
        )
    ].sort_values(
        duplicate_key
    )

    # ---------------------------------------------------------
    # Guide-level completeness table
    # ---------------------------------------------------------
    guide_hit_counts = (
        audited_hits
        .dropna(subset=["guide_id"])
        .groupby(
            ["guide_id", "gene"],
            as_index=False,
        )
        .agg(
            total_hits=(
                "guide_id",
                "size",
            ),
            genomes_with_hits=(
                "genome_normalized",
                "nunique",
            ),
            minimum_mismatches=(
                "mismatch_count_normalized",
                "min",
            ),
            maximum_mismatches=(
                "mismatch_count_normalized",
                "max",
            ),
        )
    )

    mismatch_counts = (
        audited_hits
        .dropna(subset=["guide_id"])
        .pivot_table(
            index=["guide_id", "gene"],
            columns="mismatch_count_normalized",
            values="spacer_normalized",
            aggfunc="size",
            fill_value=0,
        )
        .reset_index()
    )

    renamed_mismatch_columns = {}

    for column in mismatch_counts.columns:
        if isinstance(column, (int, float)):
            if pd.notna(column):
                renamed_mismatch_columns[column] = (
                    f"hits_mismatch_{int(column)}"
                )

    mismatch_counts = mismatch_counts.rename(
        columns=renamed_mismatch_columns
    )

    completeness = guides_clean.merge(
        guide_hit_counts,
        on=["guide_id", "gene"],
        how="left",
    )

    completeness = completeness.merge(
        mismatch_counts,
        on=["guide_id", "gene"],
        how="left",
    )

    numeric_fill_columns = [
        column
        for column in completeness.columns
        if column.startswith("hits_mismatch_")
        or column in [
            "total_hits",
            "genomes_with_hits",
        ]
    ]

    if numeric_fill_columns:
        completeness[numeric_fill_columns] = (
            completeness[numeric_fill_columns]
            .fillna(0)
            .astype(int)
        )

    completeness["cas_offinder_status"] = (
        completeness["total_hits"].apply(
            lambda value: (
                "Hits_detected"
                if value > 0
                else (
                    "Zero_hits_detected_within_"
                    "configured_search"
                )
            )
        )
    )

    completeness["explicitly_accounted_for"] = True

    completeness = completeness.sort_values(
        [
            "gene",
            "original_rank",
            "guide_id",
        ],
        na_position="last",
    )

    zero_hit_guides = completeness[
        completeness["total_hits"] == 0
    ].copy()

    # ---------------------------------------------------------
    # Other validation checks
    # ---------------------------------------------------------
    observed_genes = set(
        guides_clean["gene"].dropna()
    )

    missing_expected_genes = (
        EXPECTED_GENES - observed_genes
    )

    unexpected_genes = (
        observed_genes - EXPECTED_GENES
    )

    observed_genomes = set(
        audited_hits[
            "genome_normalized"
        ].dropna()
    )

    missing_expected_genomes = (
        EXPECTED_GENOMES - observed_genomes
    )

    unexpected_genomes = (
        observed_genomes - EXPECTED_GENOMES
    )

    invalid_mismatches = audited_hits[
        audited_hits[
            "mismatch_count_normalized"
        ].isna()
        | (
            audited_hits[
                "mismatch_count_normalized"
            ]
            < 0
        )
        | (
            audited_hits[
                "mismatch_count_normalized"
            ]
            > MAX_ALLOWED_MISMATCHES
        )
    ].copy()

    invalid_positions = audited_hits[
        audited_hits[
            "position_normalized"
        ].isna()
        | (
            audited_hits[
                "position_normalized"
            ]
            < 0
        )
    ].copy()

    invalid_strands = audited_hits[
        ~audited_hits[
            "strand_normalized"
        ].isin(
            ["+", "-"]
        )
    ].copy()

    invalid_spacer_lengths = guides_clean[
        guides_clean[
            "spacer"
        ].str.len()
        != 20
    ].copy()

    invalid_offtarget_lengths = audited_hits[
        ~audited_hits[
            "offtarget_sequence_normalized"
        ].str.len().isin(
            [20, 23]
        )
    ].copy()

    # ---------------------------------------------------------
    # Existing guide-level summary coverage
    # ---------------------------------------------------------
    summary_coverage_message = (
        "Existing guide-level summary file was not found."
    )

    existing_summary_guide_count = None
    missing_from_existing_summary: list[str] = []

    if EXISTING_GUIDE_SUMMARY_FILE.exists():
        existing_summary = pd.read_csv(
            EXISTING_GUIDE_SUMMARY_FILE
        )

        summary_guide_id_col = find_column(
            existing_summary,
            [
                "final_guide_id",
                "guide_id",
                "riskaware_guide_id",
                "selected_guide_id",
                "candidate_id",
            ],
        )

        summary_guide_ids = set(
            clean_text(
                existing_summary[
                    summary_guide_id_col
                ]
            )
        )

        selected_guide_ids = set(
            guides_clean["guide_id"]
        )

        missing_from_existing_summary = sorted(
            selected_guide_ids
            - summary_guide_ids
        )

        existing_summary_guide_count = len(
            summary_guide_ids
        )

        if missing_from_existing_summary:
            summary_coverage_message = (
                "Existing guide-level summary omits "
                "selected guides: "
                + ", ".join(
                    missing_from_existing_summary
                )
            )
        else:
            summary_coverage_message = (
                "Existing guide-level summary includes "
                "all selected guides."
            )

    # ---------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------
    completeness_file = (
        OUTPUT_DIR
        / "cas_offinder_guide_completeness_audit.csv"
    )

    zero_hit_file = (
        OUTPUT_DIR
        / "cas_offinder_zero_hit_guides.csv"
    )

    duplicate_file = (
        OUTPUT_DIR
        / "cas_offinder_duplicate_hits.csv"
    )

    unmapped_file = (
        OUTPUT_DIR
        / "cas_offinder_unmapped_hits.csv"
    )

    audited_hits_file = (
        OUTPUT_DIR
        / "cas_offinder_all_hits_audited.csv"
    )

    summary_file = (
        OUTPUT_DIR
        / "cas_offinder_qc_summary.txt"
    )

    completeness.to_csv(
        completeness_file,
        index=False,
    )

    zero_hit_guides.to_csv(
        zero_hit_file,
        index=False,
    )

    duplicate_hits.to_csv(
        duplicate_file,
        index=False,
    )

    unmapped_hits.to_csv(
        unmapped_file,
        index=False,
    )

    audited_hits.to_csv(
        audited_hits_file,
        index=False,
    )

    mismatch_distribution = (
        audited_hits[
            "mismatch_count_normalized"
        ]
        .value_counts(dropna=False)
        .sort_index()
    )

    genome_distribution = (
        audited_hits[
            "genome_normalized"
        ]
        .value_counts(dropna=False)
        .sort_index()
    )

    gene_distribution = (
        audited_hits[
            "gene"
        ]
        .value_counts(dropna=False)
        .sort_index()
    )

    strand_distribution = (
        audited_hits[
            "strand_normalized"
        ]
        .value_counts(dropna=False)
        .sort_index()
    )

    guide_count_ok = (
        len(guides_clean)
        == EXPECTED_GUIDE_COUNT
    )

    all_hits_mapped = (
        len(unmapped_hits)
        == 0
    )

    no_conflicting_ids = (
        len(conflicting_guide_ids)
        == 0
    )

    no_conflicting_genes = (
        len(conflicting_genes)
        == 0
    )

    mismatches_valid = (
        len(invalid_mismatches)
        == 0
    )

    positions_valid = (
        len(invalid_positions)
        == 0
    )

    strands_valid = (
        len(invalid_strands)
        == 0
    )

    guide_lengths_valid = (
        len(invalid_spacer_lengths)
        == 0
    )

    offtarget_lengths_valid = (
        len(invalid_offtarget_lengths)
        == 0
    )

    genomes_complete = (
        len(missing_expected_genomes)
        == 0
    )

    genes_complete = (
        len(missing_expected_genes)
        == 0
    )

    current_qc_pass = all(
        [
            guide_count_ok,
            all_hits_mapped,
            no_conflicting_ids,
            no_conflicting_genes,
            mismatches_valid,
            positions_valid,
            strands_valid,
            guide_lengths_valid,
            offtarget_lengths_valid,
            genomes_complete,
            genes_complete,
        ]
    )

    summary_lines = [
        "Cas-OFFinder Completeness Audit",
        "=" * 72,
        "",
        f"Selected guide file: {SELECTED_GUIDES_FILE}",
        f"Parsed hit file: {PARSED_HITS_FILE}",
        "",
        "Core counts",
        "-" * 72,
        (
            "Expected selected guides: "
            f"{EXPECTED_GUIDE_COUNT}"
        ),
        (
            "Observed unique selected guides: "
            f"{len(guides_clean)}"
        ),
        (
            "Total parsed hits: "
            f"{len(audited_hits)}"
        ),
        (
            "Guides with one or more hits: "
            f"{int((completeness['total_hits'] > 0).sum())}"
        ),
        (
            "Guides with zero hits: "
            f"{int((completeness['total_hits'] == 0).sum())}"
        ),
        (
            "Unmapped parsed hits: "
            f"{len(unmapped_hits)}"
        ),
        (
            "Duplicate parsed hit rows: "
            f"{len(duplicate_hits)}"
        ),
        "",
        "Expected-guide coverage",
        "-" * 72,
        (
            "Observed genes: "
            f"{sorted(observed_genes)}"
        ),
        (
            "Missing expected genes: "
            f"{sorted(missing_expected_genes)}"
        ),
        (
            "Unexpected genes: "
            f"{sorted(unexpected_genes)}"
        ),
        "",
        "Genome coverage",
        "-" * 72,
        (
            "Observed genomes: "
            f"{sorted(observed_genomes)}"
        ),
        (
            "Missing expected genomes: "
            f"{sorted(missing_expected_genomes)}"
        ),
        (
            "Unexpected genomes: "
            f"{sorted(unexpected_genomes)}"
        ),
        "",
        "Validation findings",
        "-" * 72,
        (
            "Conflicting guide IDs: "
            f"{len(conflicting_guide_ids)}"
        ),
        (
            "Conflicting gene labels: "
            f"{len(conflicting_genes)}"
        ),
        (
            "Invalid mismatch rows: "
            f"{len(invalid_mismatches)}"
        ),
        (
            "Invalid position rows: "
            f"{len(invalid_positions)}"
        ),
        (
            "Invalid strand rows: "
            f"{len(invalid_strands)}"
        ),
        (
            "Invalid selected spacer lengths: "
            f"{len(invalid_spacer_lengths)}"
        ),
        (
            "Unexpected off-target sequence lengths: "
            f"{len(invalid_offtarget_lengths)}"
        ),
        "",
        "Existing guide-level summary coverage",
        "-" * 72,
        (
            "Existing summary unique guides: "
            f"{existing_summary_guide_count}"
        ),
        summary_coverage_message,
        "",
        "Zero-hit guides",
        "-" * 72,
    ]

    if zero_hit_guides.empty:
        summary_lines.append("None")
    else:
        for row in zero_hit_guides.itertuples(
            index=False
        ):
            summary_lines.append(
                f"{row.guide_id}\t"
                f"{row.gene}\t"
                f"{row.spacer}"
            )

    summary_lines.extend(
        [
            "",
            "Hits by genome",
            "-" * 72,
            genome_distribution.to_string(),
            "",
            "Hits by target gene",
            "-" * 72,
            gene_distribution.to_string(),
            "",
            "Hits by mismatch count",
            "-" * 72,
            mismatch_distribution.to_string(),
            "",
            "Hits by strand",
            "-" * 72,
            strand_distribution.to_string(),
            "",
            "Current first-pass QC result",
            "-" * 72,
            (
                "PASS"
                if current_qc_pass
                else "REVIEW_REQUIRED"
            ),
            "",
            "Important interpretation",
            "-" * 72,
            (
                "A zero-hit guide means that no candidate site "
                "was detected in the CURRENT genome panel at "
                f"<= {MAX_ALLOWED_MISMATCHES} mismatches under "
                "the configured SpCas9 PAM pattern. It does not "
                "prove universal absence of off-targets across "
                "all strains, plasmids, or mobile genetic elements."
            ),
        ]
    )

    summary_file.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # Terminal report
    # ---------------------------------------------------------
    print("\n" + "=" * 72)
    print("AUDIT RESULT")
    print("=" * 72)

    print(
        f"Selected unique guides: "
        f"{len(guides_clean)}"
    )

    print(
        f"Parsed hits: "
        f"{len(audited_hits)}"
    )

    print(
        "Guides with hits: "
        f"{int((completeness['total_hits'] > 0).sum())}"
    )

    print(
        "Guides with zero hits: "
        f"{int((completeness['total_hits'] == 0).sum())}"
    )

    print(
        f"Unmapped hits: "
        f"{len(unmapped_hits)}"
    )

    print(
        f"Duplicate hit rows: "
        f"{len(duplicate_hits)}"
    )

    print(
        f"Invalid mismatch rows: "
        f"{len(invalid_mismatches)}"
    )

    print("\nZero-hit guides:")

    if zero_hit_guides.empty:
        print("None")
    else:
        print(
            zero_hit_guides[
                [
                    "guide_id",
                    "gene",
                    "spacer",
                    "original_rank",
                ]
            ].to_string(
                index=False
            )
        )

    print("\nHits by genome:")
    print(
        genome_distribution.to_string()
    )

    print("\nHits by mismatch count:")
    print(
        mismatch_distribution.to_string()
    )

    print("\nExisting guide-level summary:")
    print(
        summary_coverage_message
    )

    print("\nFiles written:")
    print(
        "- "
        f"{completeness_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        "- "
        f"{zero_hit_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        "- "
        f"{duplicate_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        "- "
        f"{unmapped_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        "- "
        f"{audited_hits_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        "- "
        f"{summary_file.relative_to(PROJECT_ROOT)}"
    )

    print(
        "\nCurrent first-pass QC status: "
        + (
            "PASS"
            if current_qc_pass
            else "REVIEW_REQUIRED"
        )
    )


if __name__ == "__main__":
    main()