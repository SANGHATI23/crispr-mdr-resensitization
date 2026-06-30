#!/usr/bin/env python3
"""
Create a complete guide-level Cas-OFFinder/FOTR summary.

The existing guide-level summary was generated only from guides having at
least one Cas-OFFinder hit. Consequently, four valid zero-hit guides were
omitted.

This script:

1. Loads the 20-guide Cas-OFFinder completeness audit.
2. Loads the existing 16-guide functional off-target summary.
3. Retains all existing FOTR measurements for hit-containing guides.
4. Adds the four zero-hit guides.
5. Assigns zero off-target burden only within the tested genome panel and
   configured <=4 mismatch search.
6. Produces a complete 20-guide table.
7. Generates an updated top-two-per-gene candidate table.
8. Writes a validation summary.

Inputs
------
results_cas_offinder/qc/
    cas_offinder_guide_completeness_audit.csv

results_cas_offinder/guide_level_fotr/
    guide_level_functional_offtarget_risk_summary.csv

results_final_selection/
    top5_per_gene_risk_aware_guides.csv

Outputs
-------
results_cas_offinder/guide_level_fotr_complete/
    guide_level_functional_offtarget_risk_summary_all20.csv
    zero_hit_guides_added_to_fotr.csv
    final_lowest_risk_guides_for_mic_planning_all20.csv
    complete_guide_level_fotr_summary.txt
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

COMPLETENESS_FILE = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "qc"
    / "cas_offinder_guide_completeness_audit.csv"
)

EXISTING_FOTR_FILE = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "guide_level_fotr"
    / "guide_level_functional_offtarget_risk_summary.csv"
)

SELECTED_GUIDES_FILE = (
    PROJECT_ROOT
    / "results_final_selection"
    / "top5_per_gene_risk_aware_guides.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results_cas_offinder"
    / "guide_level_fotr_complete"
)

EXPECTED_GUIDES = 20
EXPECTED_GENES = ["blaKPC", "blaNDM1", "mcr1", "mecA"]


def normalized_name(value: str) -> str:
    """Normalize a dataframe column name."""
    return (
        str(value)
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
    """Find a dataframe column from possible candidate names."""
    lookup = {
        normalized_name(column): column
        for column in dataframe.columns
    }

    for candidate in candidates:
        candidate_normalized = normalized_name(candidate)

        if candidate_normalized in lookup:
            return lookup[candidate_normalized]

    if required:
        raise KeyError(
            "Could not locate any of these columns:\n"
            f"{candidates}\n\n"
            "Available columns:\n"
            f"{list(dataframe.columns)}"
        )

    return None


def clean_text(series: pd.Series) -> pd.Series:
    """Normalize text values."""
    return series.astype(str).str.strip()


def clean_sequence(series: pd.Series) -> pd.Series:
    """Normalize DNA sequences."""
    return (
        series.astype(str)
        .str.upper()
        .str.strip()
        .str.replace(r"\s+", "", regex=True)
    )


def fail(message: str) -> None:
    """Print an error and exit."""
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    print("=" * 76)
    print("Creating complete 20-guide Cas-OFFinder/FOTR summary")
    print("=" * 76)
    print(f"Project root: {PROJECT_ROOT}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for required_file in [
        COMPLETENESS_FILE,
        EXISTING_FOTR_FILE,
        SELECTED_GUIDES_FILE,
    ]:
        if not required_file.exists():
            fail(f"Required file not found:\n{required_file}")

    completeness = pd.read_csv(COMPLETENESS_FILE)
    existing_fotr = pd.read_csv(EXISTING_FOTR_FILE)
    selected = pd.read_csv(SELECTED_GUIDES_FILE)

    print(f"\nCompleteness audit: {completeness.shape}")
    print(f"Existing FOTR summary: {existing_fotr.shape}")
    print(f"Selected guide table: {selected.shape}")

    print("\nExisting FOTR columns:")
    print(list(existing_fotr.columns))

    # ---------------------------------------------------------
    # Locate completeness-audit columns
    # ---------------------------------------------------------
    completeness_id_col = find_column(
        completeness,
        ["guide_id", "final_guide_id"],
    )

    completeness_gene_col = find_column(
        completeness,
        ["gene", "mlcb_gene"],
    )

    completeness_spacer_col = find_column(
        completeness,
        ["spacer", "mlcb_spacer", "query_spacer"],
    )

    completeness_rank_col = find_column(
        completeness,
        ["original_rank", "final_gene_rank", "guide_rank"],
    )

    completeness_total_hits_col = find_column(
        completeness,
        ["total_hits"],
    )

    completeness_genomes_col = find_column(
        completeness,
        ["genomes_with_hits"],
    )

    completeness_min_mismatch_col = find_column(
        completeness,
        ["minimum_mismatches"],
        required=False,
    )

    completeness_max_mismatch_col = find_column(
        completeness,
        ["maximum_mismatches"],
        required=False,
    )

    completeness_status_col = find_column(
        completeness,
        ["cas_offinder_status"],
    )

    # ---------------------------------------------------------
    # Locate selected-guide columns
    # ---------------------------------------------------------
    selected_id_col = find_column(
        selected,
        ["final_guide_id", "guide_id"],
    )

    selected_gene_col = find_column(
        selected,
        ["mlcb_gene", "gene"],
    )

    selected_spacer_col = find_column(
        selected,
        ["mlcb_spacer", "spacer"],
    )

    selected_rank_col = find_column(
        selected,
        ["final_gene_rank", "original_rank", "guide_rank"],
    )

    selected_base_score_col = find_column(
        selected,
        ["mlcb_base_score", "base_score"],
        required=False,
    )

    selected_annotation_risk_col = find_column(
        selected,
        [
            "annotation_derived_functional_risk_score",
            "functional_risk_score",
        ],
        required=False,
    )

    selected_predicted_risk_col = find_column(
        selected,
        [
            "predicted_high_functional_offtarget_risk_probability",
            "predicted_functional_risk",
        ],
        required=False,
    )

    selected_penalized_score_col = find_column(
        selected,
        [
            "functional_risk_penalized_score",
            "penalized_score",
        ],
        required=False,
    )

    # ---------------------------------------------------------
    # Locate existing guide-level FOTR columns
    # ---------------------------------------------------------
    fotr_id_col = find_column(
        existing_fotr,
        ["final_guide_id", "guide_id"],
    )

    fotr_gene_col = find_column(
        existing_fotr,
        ["mlcb_gene", "gene"],
    )

    # ---------------------------------------------------------
    # Build authoritative 20-guide base table
    # ---------------------------------------------------------
    complete_base = pd.DataFrame(
        {
            "final_guide_id": clean_text(
                completeness[completeness_id_col]
            ),
            "mlcb_gene": clean_text(
                completeness[completeness_gene_col]
            ),
            "mlcb_spacer": clean_sequence(
                completeness[completeness_spacer_col]
            ),
            "original_riskaware_rank": pd.to_numeric(
                completeness[completeness_rank_col],
                errors="coerce",
            ),
            "total_offtarget_hits": pd.to_numeric(
                completeness[completeness_total_hits_col],
                errors="coerce",
            ).fillna(0).astype(int),
            "genomes_with_hits": pd.to_numeric(
                completeness[completeness_genomes_col],
                errors="coerce",
            ).fillna(0).astype(int),
            "cas_offinder_status": clean_text(
                completeness[completeness_status_col]
            ),
        }
    )

    if completeness_min_mismatch_col is not None:
        complete_base["minimum_mismatches"] = pd.to_numeric(
            completeness[completeness_min_mismatch_col],
            errors="coerce",
        )
    else:
        complete_base["minimum_mismatches"] = pd.NA

    if completeness_max_mismatch_col is not None:
        complete_base["maximum_mismatches"] = pd.to_numeric(
            completeness[completeness_max_mismatch_col],
            errors="coerce",
        )
    else:
        complete_base["maximum_mismatches"] = pd.NA

    # Add selected-guide scoring fields.
    selected_normalized = pd.DataFrame(
        {
            "final_guide_id": clean_text(selected[selected_id_col]),
            "selected_gene_check": clean_text(selected[selected_gene_col]),
            "selected_spacer_check": clean_sequence(
                selected[selected_spacer_col]
            ),
            "selected_rank_check": pd.to_numeric(
                selected[selected_rank_col],
                errors="coerce",
            ),
        }
    )

    if selected_base_score_col is not None:
        selected_normalized["mlcb_base_score"] = pd.to_numeric(
            selected[selected_base_score_col],
            errors="coerce",
        )

    if selected_annotation_risk_col is not None:
        selected_normalized[
            "annotation_derived_functional_risk_score"
        ] = pd.to_numeric(
            selected[selected_annotation_risk_col],
            errors="coerce",
        )

    if selected_predicted_risk_col is not None:
        selected_normalized[
            "predicted_high_functional_offtarget_risk_probability"
        ] = pd.to_numeric(
            selected[selected_predicted_risk_col],
            errors="coerce",
        )

    if selected_penalized_score_col is not None:
        selected_normalized[
            "functional_risk_penalized_score"
        ] = pd.to_numeric(
            selected[selected_penalized_score_col],
            errors="coerce",
        )

    complete_base = complete_base.merge(
        selected_normalized,
        on="final_guide_id",
        how="left",
        validate="one_to_one",
    )

    gene_conflicts = complete_base[
        complete_base["mlcb_gene"]
        != complete_base["selected_gene_check"]
    ]

    spacer_conflicts = complete_base[
        complete_base["mlcb_spacer"]
        != complete_base["selected_spacer_check"]
    ]

    if not gene_conflicts.empty:
        fail(
            "Gene conflicts were detected between the completeness "
            "audit and selected-guide table."
        )

    if not spacer_conflicts.empty:
        fail(
            "Spacer conflicts were detected between the completeness "
            "audit and selected-guide table."
        )

    complete_base = complete_base.drop(
        columns=[
            "selected_gene_check",
            "selected_spacer_check",
            "selected_rank_check",
        ],
        errors="ignore",
    )

    # ---------------------------------------------------------
    # Normalize existing 16-guide FOTR summary
    # ---------------------------------------------------------
    fotr_normalized = existing_fotr.copy()

    fotr_normalized = fotr_normalized.rename(
        columns={
            fotr_id_col: "final_guide_id",
            fotr_gene_col: "existing_fotr_gene",
        }
    )

    fotr_normalized["final_guide_id"] = clean_text(
        fotr_normalized["final_guide_id"]
    )

    fotr_normalized["existing_fotr_gene"] = clean_text(
        fotr_normalized["existing_fotr_gene"]
    )

    # Prevent duplicate versions of fields already taken from the audit.
    fields_owned_by_base = {
        "mlcb_gene",
        "mlcb_spacer",
        "guide_id",
        "gene",
        "spacer",
        "original_rank",
        "final_gene_rank",
        "total_hits",
        "total_offtarget_hits",
        "genomes_with_hits",
        "minimum_mismatches",
        "maximum_mismatches",
        "cas_offinder_status",
    }

    columns_to_remove = []

    for column in fotr_normalized.columns:
        if (
            normalized_name(column) in fields_owned_by_base
            and column
            not in {"final_guide_id", "existing_fotr_gene"}
        ):
            columns_to_remove.append(column)

    fotr_normalized = fotr_normalized.drop(
        columns=columns_to_remove,
        errors="ignore",
    )

    # ---------------------------------------------------------
    # Merge all 20 guides with existing FOTR results
    # ---------------------------------------------------------
    complete_summary = complete_base.merge(
        fotr_normalized,
        on="final_guide_id",
        how="left",
        validate="one_to_one",
        indicator="fotr_merge_status",
    )

    existing_gene_conflicts = complete_summary[
        complete_summary["existing_fotr_gene"].notna()
        & (
            complete_summary["mlcb_gene"]
            != complete_summary["existing_fotr_gene"]
        )
    ]

    if not existing_gene_conflicts.empty:
        fail(
            "Gene conflicts were detected between the existing FOTR "
            "summary and the authoritative 20-guide table."
        )

    complete_summary = complete_summary.drop(
        columns=["existing_fotr_gene"],
        errors="ignore",
    )

    complete_summary["zero_hit_in_current_panel"] = (
        complete_summary["total_offtarget_hits"] == 0
    )

    # ---------------------------------------------------------
    # Detect important FOTR columns dynamically
    # ---------------------------------------------------------
    burden_col = find_column(
        complete_summary,
        [
            "functional_offtarget_burden_score",
            "functional_off_target_burden_score",
            "functional_burden_score",
            "burden_score",
            "fotr_burden_score",
        ],
        required=False,
    )

    high_interest_col = find_column(
        complete_summary,
        [
            "high_interest_hits",
            "high_interest_offtarget_hits",
            "high_interest_off_target_hits",
            "total_high_interest_hits",
        ],
        required=False,
    )

    recommendation_col = find_column(
        complete_summary,
        [
            "guide_level_recommendation",
            "recommendation",
            "fotr_recommendation",
            "functional_offtarget_recommendation",
        ],
        required=False,
    )

    final_rank_col = find_column(
        complete_summary,
        [
            "fotr_rank",
            "final_fotr_rank",
            "guide_level_fotr_rank",
            "functional_offtarget_rank",
        ],
        required=False,
    )

    # ---------------------------------------------------------
    # Fill zero-hit guide measurements
    # ---------------------------------------------------------
    zero_mask = complete_summary["zero_hit_in_current_panel"]

    # Only zero-hit rows receive zeros. Existing hit-containing guide
    # measurements are never overwritten.
    numeric_zero_keywords = [
        "hit",
        "burden",
        "amr",
        "resistance",
        "virulence",
        "essential",
        "mobile",
        "coding",
        "intergenic",
        "noncoding",
    ]

    protected_numeric_columns = {
        "original_riskaware_rank",
        "mlcb_base_score",
        "annotation_derived_functional_risk_score",
        "predicted_high_functional_offtarget_risk_probability",
        "functional_risk_penalized_score",
        "minimum_mismatches",
        "maximum_mismatches",
    }

    for column in complete_summary.columns:
        if column in protected_numeric_columns:
            continue

        if any(
            keyword in normalized_name(column)
            for keyword in numeric_zero_keywords
        ):
            numeric_values = pd.to_numeric(
                complete_summary[column],
                errors="coerce",
            )

            if numeric_values.notna().any() or complete_summary[
                column
            ].isna().all():
                complete_summary[column] = numeric_values
                complete_summary.loc[
                    zero_mask,
                    column,
                ] = complete_summary.loc[
                    zero_mask,
                    column,
                ].fillna(0)

    if burden_col is None:
        burden_col = "functional_offtarget_burden_score"
        complete_summary[burden_col] = pd.NA

    complete_summary[burden_col] = pd.to_numeric(
        complete_summary[burden_col],
        errors="coerce",
    )

    complete_summary.loc[
        zero_mask,
        burden_col,
    ] = 0.0

    if high_interest_col is None:
        high_interest_col = "high_interest_hits"
        complete_summary[high_interest_col] = pd.NA

    complete_summary[high_interest_col] = pd.to_numeric(
        complete_summary[high_interest_col],
        errors="coerce",
    )

    complete_summary.loc[
        zero_mask,
        high_interest_col,
    ] = 0

    if recommendation_col is None:
        recommendation_col = "guide_level_recommendation"
        complete_summary[recommendation_col] = pd.NA

    complete_summary.loc[
        zero_mask,
        recommendation_col,
    ] = (
        "Strong_candidate_zero_detected_hits_in_current_panel"
    )

    complete_summary["zero_hit_interpretation"] = ""

    complete_summary.loc[
        zero_mask,
        "zero_hit_interpretation",
    ] = (
        "No candidate off-target was detected in the current "
        "three-reference-genome panel using the configured SpCas9 "
        "PAM pattern and a maximum of four mismatches; this is not "
        "evidence of universal off-target absence."
    )

    complete_summary.loc[
        ~zero_mask,
        "zero_hit_interpretation",
    ] = "Not_applicable_hits_were_detected"

    complete_summary["fotr_record_source"] = (
        complete_summary["fotr_merge_status"].map(
            {
                "both": "Existing_hit_derived_FOTR_summary",
                "left_only": "Added_from_zero_hit_completeness_audit",
                "right_only": "Unexpected_right_only_record",
            }
        )
    )

    # ---------------------------------------------------------
    # Re-rank all guides within each gene
    # ---------------------------------------------------------
    penalized_score_col = find_column(
        complete_summary,
        ["functional_risk_penalized_score"],
        required=False,
    )

    sort_columns = [
        "mlcb_gene",
        burden_col,
        high_interest_col,
        "total_offtarget_hits",
    ]

    ascending = [
        True,
        True,
        True,
        True,
    ]

    if penalized_score_col is not None:
        sort_columns.append(penalized_score_col)
        ascending.append(False)

    sort_columns.append("original_riskaware_rank")
    ascending.append(True)

    complete_summary = complete_summary.sort_values(
        sort_columns,
        ascending=ascending,
        na_position="last",
    ).copy()

    complete_summary["complete_panel_fotr_rank"] = (
        complete_summary.groupby("mlcb_gene").cumcount() + 1
    )

    complete_summary["rank_shift_after_complete_fotr"] = (
        complete_summary["original_riskaware_rank"]
        - complete_summary["complete_panel_fotr_rank"]
    )

    # Positive rank shift means the guide moved upward.
    complete_summary["rank_shift_direction"] = (
        complete_summary[
            "rank_shift_after_complete_fotr"
        ].apply(
            lambda value: (
                "Moved_up"
                if value > 0
                else "Moved_down"
                if value < 0
                else "Unchanged"
            )
        )
    )

    # Remove the pandas merge indicator before saving.
    complete_summary = complete_summary.drop(
        columns=["fotr_merge_status"],
        errors="ignore",
    )

    # ---------------------------------------------------------
    # Validate final table
    # ---------------------------------------------------------
    unique_guides = complete_summary[
        "final_guide_id"
    ].nunique()

    duplicate_ids = complete_summary[
        complete_summary.duplicated(
            subset=["final_guide_id"],
            keep=False,
        )
    ]

    observed_genes = sorted(
        complete_summary["mlcb_gene"].dropna().unique()
    )

    per_gene_counts = (
        complete_summary["mlcb_gene"]
        .value_counts()
        .sort_index()
    )

    zero_hit_rows = complete_summary[
        complete_summary["zero_hit_in_current_panel"]
    ].copy()

    if unique_guides != EXPECTED_GUIDES:
        fail(
            f"Expected {EXPECTED_GUIDES} unique guides, "
            f"but obtained {unique_guides}."
        )

    if not duplicate_ids.empty:
        fail("Duplicate guide IDs exist in the complete summary.")

    if set(observed_genes) != set(EXPECTED_GENES):
        fail(
            "Observed genes do not match the expected genes.\n"
            f"Observed: {observed_genes}\n"
            f"Expected: {EXPECTED_GENES}"
        )

    if not all(
        per_gene_counts.get(gene, 0) == 5
        for gene in EXPECTED_GENES
    ):
        fail(
            "The completed table does not contain exactly "
            "five guides per gene."
        )

    # ---------------------------------------------------------
    # Create final top-two-per-gene table
    # ---------------------------------------------------------
    top_two = (
        complete_summary[
            complete_summary["complete_panel_fotr_rank"] <= 2
        ]
        .sort_values(
            ["mlcb_gene", "complete_panel_fotr_rank"]
        )
        .copy()
    )

    # ---------------------------------------------------------
    # Save files
    # ---------------------------------------------------------
    complete_output_file = (
        OUTPUT_DIR
        / "guide_level_functional_offtarget_risk_summary_all20.csv"
    )

    zero_hit_output_file = (
        OUTPUT_DIR
        / "zero_hit_guides_added_to_fotr.csv"
    )

    top_two_output_file = (
        OUTPUT_DIR
        / "final_lowest_risk_guides_for_mic_planning_all20.csv"
    )

    summary_output_file = (
        OUTPUT_DIR
        / "complete_guide_level_fotr_summary.txt"
    )

    complete_summary.to_csv(
        complete_output_file,
        index=False,
    )

    zero_hit_rows.to_csv(
        zero_hit_output_file,
        index=False,
    )

    top_two.to_csv(
        top_two_output_file,
        index=False,
    )

    summary_lines = [
        "Complete Guide-Level Cas-OFFinder/FOTR Summary",
        "=" * 76,
        "",
        "Input files",
        "-" * 76,
        f"Completeness audit: {COMPLETENESS_FILE}",
        f"Existing FOTR summary: {EXISTING_FOTR_FILE}",
        f"Selected guides: {SELECTED_GUIDES_FILE}",
        "",
        "Completion results",
        "-" * 76,
        f"Total guide rows: {len(complete_summary)}",
        f"Unique guide IDs: {unique_guides}",
        f"Observed genes: {observed_genes}",
        f"Zero-hit guides added: {len(zero_hit_rows)}",
        f"Guides with detected hits: {int((~zero_mask).sum())}",
        "",
        "Guide count per gene",
        "-" * 76,
        per_gene_counts.to_string(),
        "",
        "Zero-hit guides added",
        "-" * 76,
    ]

    for row in zero_hit_rows.itertuples(index=False):
        summary_lines.append(
            f"{row.final_guide_id}\t"
            f"{row.mlcb_gene}\t"
            f"{row.mlcb_spacer}"
        )

    summary_lines.extend(
        [
            "",
            "Top two complete-panel FOTR guides per gene",
            "-" * 76,
        ]
    )

    for row in top_two.itertuples(index=False):
        summary_lines.append(
            f"{row.mlcb_gene}\t"
            f"rank={row.complete_panel_fotr_rank}\t"
            f"{row.final_guide_id}\t"
            f"hits={row.total_offtarget_hits}\t"
            f"burden={getattr(row, burden_col)}"
        )

    summary_lines.extend(
        [
            "",
            "Validation",
            "-" * 76,
            "Expected total guides: 20",
            f"Observed total guides: {unique_guides}",
            "Expected guides per gene: 5",
            (
                "All genes contain five guides: "
                + str(
                    all(
                        per_gene_counts.get(gene, 0) == 5
                        for gene in EXPECTED_GENES
                    )
                )
            ),
            f"Duplicate guide IDs: {len(duplicate_ids)}",
            "Status: PASS",
            "",
            "Scientific interpretation",
            "-" * 76,
            (
                "The four added guides had no candidate off-target "
                "hits in the current three-reference-genome panel "
                "under the configured SpCas9 NGG search with a "
                "maximum of four mismatches. They must not be "
                "described as universally off-target-free. Their "
                "zero-hit status must be re-evaluated after expanding "
                "the analysis to clinical strains, plasmids, SCCmec, "
                "integrons, resistance islands, and other mobile "
                "genetic elements."
            ),
        ]
    )

    summary_output_file.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    # ---------------------------------------------------------
    # Terminal output
    # ---------------------------------------------------------
    print("\n" + "=" * 76)
    print("COMPLETE GUIDE-LEVEL SUMMARY RESULT")
    print("=" * 76)

    print(f"Total guide rows: {len(complete_summary)}")
    print(f"Unique guides: {unique_guides}")
    print(f"Zero-hit guides added: {len(zero_hit_rows)}")
    print(f"Duplicate guide IDs: {len(duplicate_ids)}")

    print("\nGuide count per gene:")
    print(per_gene_counts.to_string())

    print("\nZero-hit guides added:")
    print(
        zero_hit_rows[
            [
                "final_guide_id",
                "mlcb_gene",
                "mlcb_spacer",
                "total_offtarget_hits",
                burden_col,
                recommendation_col,
            ]
        ].to_string(index=False)
    )

    print("\nUpdated top two guides per gene:")
    display_columns = [
        "mlcb_gene",
        "complete_panel_fotr_rank",
        "final_guide_id",
        "mlcb_spacer",
        "total_offtarget_hits",
        burden_col,
        high_interest_col,
        recommendation_col,
    ]

    print(
        top_two[display_columns].to_string(index=False)
    )

    print("\nFiles written:")
    print(
        f"- {complete_output_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        f"- {zero_hit_output_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        f"- {top_two_output_file.relative_to(PROJECT_ROOT)}"
    )
    print(
        f"- {summary_output_file.relative_to(PROJECT_ROOT)}"
    )

    print("\nComplete 20-guide FOTR summary status: PASS")


if __name__ == "__main__":
    main()