#!/usr/bin/env python3
"""
01_prepare_mlcb_dataset.py

Purpose
-------
Prepare a clean MLCB-ready guide-level dataset for the FOTR-CRISPR + DNABERT-2 paper.

This script:
1. Selects the full 145-guide comparison table as the main MLCB dataset.
2. Standardizes guide/gene/PAM/score columns.
3. Creates DNABERT-ready sequence fields.
4. Creates weak labels for MLCB experiments.
5. Writes clean output files to results_mlcb/tables/.

Main expected input
-------------------
results/final_comparison_table.csv

Expected output
---------------
results_mlcb/tables/mlcb_guide_dataset.csv
results_mlcb/tables/mlcb_guide_dataset_summary.csv
results_mlcb/tables/mlcb_dataset_preparation_log.txt
"""

from pathlib import Path
import pandas as pd
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_DIR = PROJECT_ROOT / "results_mlcb" / "tables"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_DATASET = OUTPUT_DIR / "mlcb_guide_dataset.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "mlcb_guide_dataset_summary.csv"
OUTPUT_LOG = OUTPUT_DIR / "mlcb_dataset_preparation_log.txt"


# IMPORTANT:
# Keep the full 145-guide table first.
# Do not put final_tables/table3_top_fotr_v2_guides_by_gene.csv first,
# because that file has only 4 final rows.
CANDIDATE_INPUT_FILES = [
    PROJECT_ROOT / "results" / "final_comparison_table.csv",
    PROJECT_ROOT / "results" / "cross_model_comparison.csv",
    PROJECT_ROOT / "results_case_study" / "meca_case_study_all_guides.csv",
    PROJECT_ROOT / "final_tables" / "table3_top_fotr_v2_guides_by_gene.csv",
    PROJECT_ROOT / "results_external" / "crisot" / "merged" / "crisot_guide_level_summary.csv",
    PROJECT_ROOT / "results_case_study" / "meca_case_study_top_candidates.csv",
]


def log(message: str, lines: list[str]) -> None:
    print(message)
    lines.append(message)


def find_input_file(log_lines: list[str]) -> Path:
    """
    Find the best available CSV file for MLCB guide-level preparation.
    """
    log("Searching candidate input files...", log_lines)

    existing = []
    for path in CANDIDATE_INPUT_FILES:
        if path.exists():
            existing.append(path)
            log(f"FOUND: {path.relative_to(PROJECT_ROOT)}", log_lines)
        else:
            log(f"missing: {path.relative_to(PROJECT_ROOT)}", log_lines)

    if not existing:
        raise FileNotFoundError(
            "No candidate input CSV found. Please check results/, final_tables/, "
            "results_external/, or results_case_study/."
        )

    selected = existing[0]
    log(f"\nSelected input file: {selected.relative_to(PROJECT_ROOT)}", log_lines)
    return selected


def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize column names to lower_snake_case style.
    """
    df = df.copy()
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "_", regex=False)
        .str.replace("/", "_", regex=False)
        .str.replace("__", "_", regex=False)
        .str.lower()
    )
    return df


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """
    Return the first matching column from a candidate list.
    """
    for col in candidates:
        if col in df.columns:
            return col
    return None


def calculate_gc_content(seq: str) -> float:
    """
    Calculate GC percentage for a DNA sequence.
    """
    seq = str(seq).upper().strip()
    if len(seq) == 0:
        return np.nan

    gc = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


def create_standard_columns(df: pd.DataFrame, log_lines: list[str]) -> pd.DataFrame:
    """
    Create standard MLCB columns expected by later scripts.
    """
    df = df.copy()

    gene_col = first_existing_column(
        df,
        [
            "gene",
            "target_gene",
            "resistance_gene",
        ],
    )

    spacer_col = first_existing_column(
        df,
        [
            "spacer",
            "top_spacer",
            "guide",
            "guide_sequence",
            "grna",
            "sgrna",
        ],
    )

    pam_col = first_existing_column(
        df,
        [
            "pam",
            "pam_sequence",
        ],
    )

    position_col = first_existing_column(
        df,
        [
            "position_1based",
            "position",
            "target_position",
        ],
    )

    gc_col = first_existing_column(
        df,
        [
            "gc_content",
            "gc",
            "gc_percent",
        ],
    )

    region_col = first_existing_column(
        df,
        [
            "target_region_bin",
            "target_region",
            "region",
        ],
    )

    rna_col = first_existing_column(
        df,
        [
            "rna_accessibility_class",
            "rna_accessibility",
            "accessibility_class",
        ],
    )

    escape_col = first_existing_column(
        df,
        [
            "escape_risk_class",
            "crisot_escape_risk_class",
            "crisot_risk_class",
        ],
    )

    score_col = first_existing_column(
        df,
        [
            "fotr_v2_priority_score",
            "fotr_v1_priority_score",
            "final_fotr_score",
            "final_score",
            "priority_score",
            "unseen_gene_priority_score_v1",
            "clinical_pareto_support_score",
            "simple_on_target_proxy",
            "on_target_score",
            "rs3_score",
        ],
    )

    label_col = first_existing_column(
        df,
        [
            "fotr_v2_recommendation_status",
            "final_recommendation",
            "recommendation",
            "classification",
            "unseen_gene_recommendation_v1",
            "final_clinical_pareto_label",
            "fotr_recommendation",
            "recommended_for_final_table",
        ],
    )

    required_map = {
        "gene": gene_col,
        "spacer": spacer_col,
    }

    missing_required = [std for std, real in required_map.items() if real is None]
    if missing_required:
        raise ValueError(
            f"Missing required columns: {missing_required}. "
            f"Available columns: {list(df.columns)}"
        )

    df["mlcb_gene"] = df[gene_col].astype(str).str.strip()
    df["mlcb_spacer"] = df[spacer_col].astype(str).str.upper().str.strip()

    if pam_col:
        df["mlcb_pam"] = df[pam_col].astype(str).str.upper().str.strip()
    else:
        # final_comparison_table.csv may not contain PAM.
        # Use NNN placeholder so DNABERT sequence creation still works.
        df["mlcb_pam"] = "NNN"
        log("No PAM column found. mlcb_pam set to NNN placeholder.", log_lines)

    if position_col:
        df["mlcb_position_1based"] = pd.to_numeric(df[position_col], errors="coerce")
    else:
        df["mlcb_position_1based"] = np.nan
        log("No position column found. mlcb_position_1based set to NaN.", log_lines)

    if gc_col:
        df["mlcb_gc_content"] = pd.to_numeric(df[gc_col], errors="coerce")
    else:
        df["mlcb_gc_content"] = df["mlcb_spacer"].apply(calculate_gc_content)
        log("No GC column found. Calculated GC content from spacer.", log_lines)

    if region_col:
        df["mlcb_target_region_bin"] = df[region_col].astype(str)
    else:
        df["mlcb_target_region_bin"] = "unknown"
        log("No target-region column found. mlcb_target_region_bin set to unknown.", log_lines)

    if rna_col:
        df["mlcb_rna_accessibility_class"] = df[rna_col].astype(str)
    else:
        df["mlcb_rna_accessibility_class"] = "unknown"
        log("No RNA accessibility column found. mlcb_rna_accessibility_class set to unknown.", log_lines)

    if escape_col:
        df["mlcb_escape_risk_class"] = df[escape_col].astype(str)
    else:
        df["mlcb_escape_risk_class"] = "unknown"
        log("No escape-risk column found. mlcb_escape_risk_class set to unknown.", log_lines)

    if score_col:
        df["mlcb_base_score"] = pd.to_numeric(df[score_col], errors="coerce")
        log(f"Using score column: {score_col}", log_lines)
    else:
        df["mlcb_base_score"] = np.nan
        log("No score column found. mlcb_base_score set to NaN.", log_lines)

    if label_col:
        df["mlcb_original_recommendation"] = df[label_col].astype(str)
        log(f"Using recommendation/classification column: {label_col}", log_lines)
    else:
        df["mlcb_original_recommendation"] = "unknown"
        log("No recommendation column found. mlcb_original_recommendation set to unknown.", log_lines)

    df["mlcb_spacer_pam_sequence"] = df["mlcb_spacer"] + df["mlcb_pam"]

    return df


def create_dnabert_input_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create sequence fields for DNABERT-2.

    Current version:
    - spacer_only_sequence = 20 nt spacer
    - spacer_pam_sequence = spacer + PAM
    - dnabert_sequence = spacer + PAM

    Later upgrade:
    Replace dnabert_sequence with 60-100 bp guide-centered genomic context
    once FASTA extraction is added.
    """
    df = df.copy()

    df["spacer_only_sequence"] = df["mlcb_spacer"]
    df["spacer_pam_sequence"] = df["mlcb_spacer_pam_sequence"]

    df["dnabert_sequence"] = df["spacer_pam_sequence"]
    df["dnabert_sequence_type"] = "spacer_plus_pam_or_placeholder"

    return df


def create_weak_labels(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create weak labels for MLCB experiments.

    These are not experimental labels.
    They are derived labels for representation/probing analysis.
    """
    df = df.copy()

    recommendation = df["mlcb_original_recommendation"].astype(str).str.lower()

    # High-priority guide label.
    # If classification exists, we treat strong/selected/final/recommended labels as 1.
    df["label_high_priority_guide"] = np.where(
        recommendation.str.contains(
            "high|final|recommended|pareto|candidate|selected|good|top",
            regex=True,
        ),
        1,
        0,
    )

    # If score exists, also mark top quartile as high priority.
    # This helps when recommendation label is missing or generic.
    if df["mlcb_base_score"].notna().sum() > 0:
        score_threshold = df["mlcb_base_score"].quantile(0.75)
        df["label_high_priority_by_score_top_quartile"] = np.where(
            df["mlcb_base_score"] >= score_threshold,
            1,
            0,
        )
    else:
        df["label_high_priority_by_score_top_quartile"] = -1

    # RNA structural risk weak label.
    rna = df["mlcb_rna_accessibility_class"].astype(str).str.lower()
    df["label_rna_structural_risk"] = np.select(
        [
            rna.str.contains("accessible"),
            rna.str.contains("moderate"),
            rna.str.contains("risky|risk|inaccessible|structural"),
        ],
        [0, 1, 2],
        default=-1,
    )

    # Escape risk weak label.
    escape = df["mlcb_escape_risk_class"].astype(str).str.lower()
    df["label_escape_risk"] = np.select(
        [
            escape.str.contains("low"),
            escape.str.contains("moderate|medium"),
            escape.str.contains("high|risky|risk"),
        ],
        [0, 1, 2],
        default=-1,
    )

    # Placeholder functional severity label.
    # If target-region information exists:
    # early coding = high because disrupting early coding regions is often more functionally severe.
    region = df["mlcb_target_region_bin"].astype(str).str.lower()

    df["weak_functional_severity_label"] = np.select(
        [
            region.str.contains("early"),
            region.str.contains("middle|mid"),
            region.str.contains("late|unknown|intergenic"),
        ],
        [2, 1, 0],
        default=1,
    )

    df["weak_functional_severity_name"] = df["weak_functional_severity_label"].map(
        {
            0: "low",
            1: "moderate",
            2: "high",
        }
    )

    return df


def remove_bad_rows(df: pd.DataFrame, log_lines: list[str]) -> pd.DataFrame:
    """
    Keep rows with valid DNA guide sequences.
    """
    before = len(df)

    df = df.copy()

    df = df[df["mlcb_spacer"].notna()]
    df = df[df["mlcb_spacer"].str.len() >= 18]
    df = df[df["mlcb_spacer"].str.contains(r"^[ACGTN]+$", regex=True)]

    after = len(df)
    log(f"Rows before sequence cleaning: {before}", log_lines)
    log(f"Rows after sequence cleaning:  {after}", log_lines)

    return df


def create_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a compact summary table.
    """
    summary_rows = []

    summary_rows.append({"metric": "total_rows", "value": len(df)})
    summary_rows.append({"metric": "unique_genes", "value": df["mlcb_gene"].nunique()})
    summary_rows.append({"metric": "unique_spacers", "value": df["mlcb_spacer"].nunique()})

    for gene, count in df["mlcb_gene"].value_counts().items():
        summary_rows.append({"metric": f"guides_for_gene_{gene}", "value": count})

    for label, count in df["label_high_priority_guide"].value_counts(dropna=False).items():
        summary_rows.append({"metric": f"high_priority_text_label_{label}", "value": count})

    for label, count in df["label_high_priority_by_score_top_quartile"].value_counts(dropna=False).items():
        summary_rows.append({"metric": f"high_priority_score_quartile_label_{label}", "value": count})

    for label, count in df["weak_functional_severity_name"].value_counts(dropna=False).items():
        summary_rows.append({"metric": f"weak_functional_severity_{label}", "value": count})

    return pd.DataFrame(summary_rows)


def main() -> None:
    log_lines = []

    log("Starting MLCB dataset preparation...", log_lines)
    log(f"Project root: {PROJECT_ROOT}", log_lines)

    input_file = find_input_file(log_lines)

    df = pd.read_csv(input_file)
    log(f"Loaded input shape: {df.shape}", log_lines)

    df = normalize_column_names(df)
    log(f"Normalized columns: {list(df.columns)}", log_lines)

    df = create_standard_columns(df, log_lines)
    df = create_dnabert_input_sequences(df)
    df = create_weak_labels(df)
    df = remove_bad_rows(df, log_lines)

    final_columns = [
        "mlcb_gene",
        "mlcb_spacer",
        "mlcb_pam",
        "mlcb_spacer_pam_sequence",
        "spacer_only_sequence",
        "spacer_pam_sequence",
        "dnabert_sequence",
        "dnabert_sequence_type",
        "mlcb_position_1based",
        "mlcb_gc_content",
        "mlcb_target_region_bin",
        "mlcb_rna_accessibility_class",
        "mlcb_escape_risk_class",
        "mlcb_base_score",
        "mlcb_original_recommendation",
        "label_high_priority_guide",
        "label_high_priority_by_score_top_quartile",
        "label_rna_structural_risk",
        "label_escape_risk",
        "weak_functional_severity_label",
        "weak_functional_severity_name",
    ]

    existing_final_columns = [col for col in final_columns if col in df.columns]
    output_df = df[existing_final_columns].copy()

    output_df.to_csv(OUTPUT_DATASET, index=False)

    summary_df = create_summary(output_df)
    summary_df.to_csv(OUTPUT_SUMMARY, index=False)

    log(f"\nWrote dataset: {OUTPUT_DATASET.relative_to(PROJECT_ROOT)}", log_lines)
    log(f"Wrote summary: {OUTPUT_SUMMARY.relative_to(PROJECT_ROOT)}", log_lines)
    log(f"Final output shape: {output_df.shape}", log_lines)

    with open(OUTPUT_LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines))

    log(f"Wrote log: {OUTPUT_LOG.relative_to(PROJECT_ROOT)}", log_lines)

    print("\nPreview:")
    print(output_df.head(10).to_string(index=False))

    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()