#!/usr/bin/env python3
"""
Create manuscript-ready final tables for the FOTR-CRISPR BIBM paper.

This script reads the completed FOTR-CRISPR result files and generates
clean CSV tables in final_tables/.

Generated tables:
1. table1_dataset_summary.csv
2. table2_rna_accessibility_by_gene.csv
3. table3_top_fotr_v2_guides_by_gene.csv
4. table4_recommendation_distribution_by_gene.csv
5. table5_ablation_summary.csv
6. table6_meca_case_study_summary.csv

Run from repository root:

    python analysis/manuscript_outputs/create_final_tables.py
"""

from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FINAL_TABLES_DIR = ROOT / "final_tables"
FINAL_TABLES_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------

RNA_FILE = ROOT / "results_rna" / "all_guides_rna_structure.csv"
FOTR_V2_FILE = (
    ROOT
    / "results_fotr_v2"
    / "all_guides_fotr_v2_functional_context_recommended.csv"
)
ABLATION_FILE = ROOT / "results_ablation" / "fotr_v2_ablation_summary.csv"
MECA_TOP_FILE = ROOT / "results_case_study" / "meca_case_study_top_candidates.csv"
MECA_SUMMARY_FILE = ROOT / "results_case_study" / "meca_case_study_summary.txt"


def require_file(path: Path) -> None:
    """Stop clearly if an expected input file is missing."""
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")


def normalize_gene_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the gene column name across outputs.

    Different analysis scripts may use gene, Gene, gene_name, or gene_normalized.
    This function creates a consistent 'gene' column.
    """
    candidate_cols = ["gene", "Gene", "gene_name", "target_gene", "gene_normalized"]

    for col in candidate_cols:
        if col in df.columns:
            df = df.copy()
            df["gene"] = df[col]
            return df

    raise ValueError(
        f"No recognizable gene column found. Available columns: {list(df.columns)}"
    )


def find_column(df: pd.DataFrame, candidates: list[str], required: bool = True):
    """
    Find the first matching column from a list of candidates.
    """
    for col in candidates:
        if col in df.columns:
            return col

    if required:
        raise ValueError(
            f"None of the candidate columns were found: {candidates}. "
            f"Available columns: {list(df.columns)}"
        )

    return None


def save_table(df: pd.DataFrame, filename: str) -> None:
    """Save a table to final_tables/."""
    output_path = FINAL_TABLES_DIR / filename
    df.to_csv(output_path, index=False)
    print(f"Wrote: {output_path.relative_to(ROOT)}")


# ---------------------------------------------------------------------
# Table 1: Dataset summary
# ---------------------------------------------------------------------

def create_table1_dataset_summary(rna_df: pd.DataFrame, fotr_df: pd.DataFrame) -> None:
    rna_df = normalize_gene_column(rna_df)
    fotr_df = normalize_gene_column(fotr_df)

    spacer_col = find_column(
        fotr_df,
        ["spacer", "Spacer", "guide_sequence", "Guide", "target_sequence"],
        required=False,
    )

    pam_col = find_column(
        fotr_df,
        ["pam", "PAM"],
        required=False,
    )

    rows = []

    for gene, group in fotr_df.groupby("gene"):
        row = {
            "gene": gene,
            "fotr_v2_guide_rows": len(group),
        }

        if spacer_col:
            row["unique_spacers"] = group[spacer_col].nunique()

        if pam_col:
            row["unique_pams"] = group[pam_col].nunique()

        rows.append(row)

    table = pd.DataFrame(rows).sort_values("gene")

    total_row = {
        "gene": "TOTAL",
        "fotr_v2_guide_rows": len(fotr_df),
    }

    if spacer_col:
        total_row["unique_spacers"] = fotr_df[spacer_col].nunique()

    if pam_col:
        total_row["unique_pams"] = fotr_df[pam_col].nunique()

    table = pd.concat([table, pd.DataFrame([total_row])], ignore_index=True)

    save_table(table, "table1_dataset_summary.csv")


# ---------------------------------------------------------------------
# Table 2: RNA accessibility by gene
# ---------------------------------------------------------------------

def create_table2_rna_accessibility_by_gene(rna_df: pd.DataFrame) -> None:
    rna_df = normalize_gene_column(rna_df)

    rna_class_col = find_column(
        rna_df,
        [
            "rna_accessibility_class",
            "accessibility_class",
            "RNA_Class",
            "rna_class",
        ],
    )

    count_table = (
        rna_df.groupby(["gene", rna_class_col])
        .size()
        .reset_index(name="guide_count")
    )

    pivot = count_table.pivot(
        index="gene",
        columns=rna_class_col,
        values="guide_count",
    ).fillna(0)

    for col in ["Accessible", "Moderate", "Structurally_Risky"]:
        if col not in pivot.columns:
            pivot[col] = 0

    pivot = pivot[["Accessible", "Moderate", "Structurally_Risky"]]
    pivot["total_guides"] = pivot.sum(axis=1)

    for col in ["Accessible", "Moderate", "Structurally_Risky"]:
        pivot[f"{col}_percent"] = (
            pivot[col] / pivot["total_guides"] * 100
        ).round(2)

    table = pivot.reset_index().sort_values("gene")

    save_table(table, "table2_rna_accessibility_by_gene.csv")


# ---------------------------------------------------------------------
# Table 3: Top FOTR v2 guides by gene
# ---------------------------------------------------------------------

def create_table3_top_fotr_v2_guides_by_gene(fotr_df: pd.DataFrame) -> None:
    fotr_df = normalize_gene_column(fotr_df)

    score_col = find_column(
        fotr_df,
        [
            "fotr_v2_priority_score",
            "fotr_v2_score",
            "FOTR_v2_score",
            "priority_score_v2",
        ],
    )

    spacer_col = find_column(
        fotr_df,
        ["spacer", "Spacer", "guide_sequence", "Guide", "target_sequence"],
    )

    pam_col = find_column(fotr_df, ["pam", "PAM"], required=False)

    rna_class_col = find_column(
        fotr_df,
        [
            "rna_accessibility_class",
            "accessibility_class",
            "RNA_Class",
            "rna_class",
        ],
        required=False,
    )

    target_region_col = find_column(
        fotr_df,
        ["target_region_bin", "target_region", "region_bin"],
        required=False,
    )

    recommendation_col = find_column(
        fotr_df,
        [
            "fotr_v2_recommendation_status",
            "recommendation_status",
            "recommendation",
            "final_recommendation",
            "priority_label",
        ],
        required=False,
    )

    recommended_flag_col = find_column(
        fotr_df,
        [
            "recommended_for_final_table",
            "recommended",
            "include_in_final_table",
        ],
        required=False,
    )

    # Prefer final-table recommended guides when available.
    ranking_df = fotr_df.copy()

    if recommended_flag_col:
        recommended_values = ranking_df[recommended_flag_col].astype(str).str.lower()
        recommended_df = ranking_df[
            recommended_values.isin(["true", "1", "yes", "y"])
        ].copy()

        if not recommended_df.empty:
            ranking_df = recommended_df

    columns = ["gene", spacer_col, score_col]

    if pam_col:
        columns.append(pam_col)
    if rna_class_col:
        columns.append(rna_class_col)
    if target_region_col:
        columns.append(target_region_col)
    if recommendation_col:
        columns.append(recommendation_col)
    if recommended_flag_col:
        columns.append(recommended_flag_col)

    top_rows = (
        ranking_df.sort_values(["gene", score_col], ascending=[True, False])
        .groupby("gene")
        .head(1)
        .loc[:, columns]
        .copy()
    )

    rename_map = {
        spacer_col: "top_spacer",
        score_col: "fotr_v2_priority_score",
    }

    if pam_col:
        rename_map[pam_col] = "pam"
    if rna_class_col:
        rename_map[rna_class_col] = "rna_accessibility_class"
    if target_region_col:
        rename_map[target_region_col] = "target_region"
    if recommendation_col:
        rename_map[recommendation_col] = "fotr_v2_recommendation_status"
    if recommended_flag_col:
        rename_map[recommended_flag_col] = "recommended_for_final_table"

    top_rows = top_rows.rename(columns=rename_map)
    top_rows["fotr_v2_priority_score"] = top_rows["fotr_v2_priority_score"].round(2)

    save_table(top_rows, "table3_top_fotr_v2_guides_by_gene.csv")


# ---------------------------------------------------------------------
# Table 4: Recommendation distribution by gene
# ---------------------------------------------------------------------

def create_table4_recommendation_distribution_by_gene(fotr_df: pd.DataFrame) -> None:
    fotr_df = normalize_gene_column(fotr_df)

    recommendation_col = find_column(
        fotr_df,
        [
            "fotr_v2_recommendation_status",
            "recommendation_status",
            "recommendation",
            "final_recommendation",
            "priority_label",
        ],
    )

    count_table = (
        fotr_df.groupby(["gene", recommendation_col])
        .size()
        .reset_index(name="guide_count")
    )

    pivot = count_table.pivot(
        index="gene",
        columns=recommendation_col,
        values="guide_count",
    ).fillna(0)

    preferred_order = [
        "High_priority_recommended",
        "Moderate_priority_candidate",
        "Not_recommended_structural_risk",
    ]

    existing_preferred = [col for col in preferred_order if col in pivot.columns]
    remaining = [col for col in pivot.columns if col not in existing_preferred]
    pivot = pivot[existing_preferred + remaining]

    pivot["total_guides"] = pivot.sum(axis=1)

    table = pivot.reset_index().sort_values("gene")

    save_table(table, "table4_recommendation_distribution_by_gene.csv")


# ---------------------------------------------------------------------
# Table 5: Ablation summary
# ---------------------------------------------------------------------

def create_table5_ablation_summary(ablation_df: pd.DataFrame) -> None:
    ablation_df = normalize_gene_column(ablation_df)

    removed_col = find_column(
        ablation_df,
        [
            "removed_component",
            "removed_feature_group",
            "ablation",
            "ablated_component",
        ],
    )

    mean_shift_col = find_column(
        ablation_df,
        [
            "mean_absolute_rank_shift",
            "mean_abs_rank_shift",
            "mean_rank_shift",
        ],
    )

    max_shift_col = find_column(
        ablation_df,
        [
            "max_rank_shift",
            "maximum_rank_shift",
            "max_absolute_rank_shift",
        ],
        required=False,
    )

    spearman_col = find_column(
        ablation_df,
        [
            "spearman_correlation",
            "spearman_rank_correlation",
            "spearman",
        ],
        required=False,
    )

    top10_col = find_column(
        ablation_df,
        [
            "top10_jaccard",
            "top_10_jaccard",
            "top10_jaccard_similarity",
            "top_10_overlap",
            "top10_overlap",
            "top_10_overlap_count",
        ],
        required=False,
    )

    idx = ablation_df.groupby("gene")[mean_shift_col].idxmax()
    table = ablation_df.loc[idx].copy()

    keep_cols = ["gene", removed_col, mean_shift_col]

    if max_shift_col:
        keep_cols.append(max_shift_col)
    if spearman_col:
        keep_cols.append(spearman_col)
    if top10_col:
        keep_cols.append(top10_col)

    table = table[keep_cols].copy()

    rename_map = {
        removed_col: "most_influential_removed_component",
        mean_shift_col: "mean_absolute_rank_shift",
    }

    if max_shift_col:
        rename_map[max_shift_col] = "max_rank_shift"
    if spearman_col:
        rename_map[spearman_col] = "spearman_rank_correlation"
    if top10_col:
        rename_map[top10_col] = "top10_jaccard_or_overlap"

    table = table.rename(columns=rename_map)

    numeric_cols = table.select_dtypes(include="number").columns
    table[numeric_cols] = table[numeric_cols].round(4)

    table = table.sort_values("gene")

    save_table(table, "table5_ablation_summary.csv")


# ---------------------------------------------------------------------
# Table 6: mecA case study summary
# ---------------------------------------------------------------------

def create_table6_meca_case_study_summary() -> None:
    require_file(MECA_TOP_FILE)

    meca_df = pd.read_csv(MECA_TOP_FILE)
    meca_df = normalize_gene_column(meca_df)

    score_col = find_column(
        meca_df,
        [
            "fotr_v2_priority_score",
            "fotr_v2_score",
            "FOTR_v2_score",
            "priority_score_v2",
        ],
        required=False,
    )

    spacer_col = find_column(
        meca_df,
        ["spacer", "Spacer", "guide_sequence", "Guide", "target_sequence"],
        required=False,
    )

    pam_col = find_column(meca_df, ["pam", "PAM"], required=False)

    position_col = find_column(
        meca_df,
        ["position", "Position", "guide_position", "start", "start_position"],
        required=False,
    )

    rna_class_col = find_column(
        meca_df,
        [
            "rna_accessibility_class",
            "accessibility_class",
            "RNA_Class",
            "rna_class",
        ],
        required=False,
    )

    target_region_col = find_column(
        meca_df,
        ["target_region_bin", "target_region", "region_bin"],
        required=False,
    )

    recommendation_col = find_column(
        meca_df,
        [
            "fotr_v2_recommendation_status",
            "recommendation_status",
            "recommendation",
            "final_recommendation",
            "priority_label",
        ],
        required=False,
    )

    recommended_flag_col = find_column(
        meca_df,
        [
            "recommended_for_final_table",
            "recommended",
            "include_in_final_table",
        ],
        required=False,
    )

    # Prefer high-priority recommended guides if present.
    top_df = meca_df.copy()

    if recommendation_col:
        high_priority = top_df[
            top_df[recommendation_col].astype(str) == "High_priority_recommended"
        ].copy()

        if not high_priority.empty:
            top_df = high_priority

    if score_col:
        top_df = top_df.sort_values(score_col, ascending=False)

    top = top_df.head(1).copy()

    output = pd.DataFrame(
        [
            {
                "gene": top["gene"].iloc[0] if "gene" in top.columns else "mecA",
                "top_spacer": top[spacer_col].iloc[0] if spacer_col else "",
                "pam": top[pam_col].iloc[0] if pam_col else "",
                "position": top[position_col].iloc[0] if position_col else "",
                "rna_accessibility_class": top[rna_class_col].iloc[0]
                if rna_class_col
                else "",
                "target_region": top[target_region_col].iloc[0]
                if target_region_col
                else "",
                "fotr_v2_priority_score": round(float(top[score_col].iloc[0]), 2)
                if score_col
                else "",
                "fotr_v2_recommendation_status": top[recommendation_col].iloc[0]
                if recommendation_col
                else "",
                "recommended_for_final_table": top[recommended_flag_col].iloc[0]
                if recommended_flag_col
                else "",
                "case_study_interpretation": (
                    "mecA showed the strongest RNA-structure filtering effect; "
                    "most guide rows were structurally risky, while the top final "
                    "candidate remained accessible and targeted an early coding region."
                ),
            }
        ]
    )

    save_table(output, "table6_meca_case_study_summary.csv")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    require_file(RNA_FILE)
    require_file(FOTR_V2_FILE)
    require_file(ABLATION_FILE)
    require_file(MECA_TOP_FILE)

    print("Reading completed FOTR-CRISPR result files...")

    rna_df = pd.read_csv(RNA_FILE)
    fotr_df = pd.read_csv(FOTR_V2_FILE)
    ablation_df = pd.read_csv(ABLATION_FILE)

    print("Creating manuscript-ready final tables...")

    create_table1_dataset_summary(rna_df, fotr_df)
    create_table2_rna_accessibility_by_gene(rna_df)
    create_table3_top_fotr_v2_guides_by_gene(fotr_df)
    create_table4_recommendation_distribution_by_gene(fotr_df)
    create_table5_ablation_summary(ablation_df)
    create_table6_meca_case_study_summary()

    print("\nDone. Final tables are available in:")
    print(f"  {FINAL_TABLES_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()