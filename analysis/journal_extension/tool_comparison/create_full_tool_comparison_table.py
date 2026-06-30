from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[3]

INPUT_TOOL_TABLE = ROOT / "results" / "final_comparison_table.csv"
INPUT_FOTR_V2 = ROOT / "results_fotr_v2" / "all_guides_fotr_v2_functional_context_recommended.csv"

OUT_DIR = ROOT / "results_journal_extension" / "tool_comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "full_tool_comparison_table.csv"
OUT_SUMMARY = OUT_DIR / "tool_comparison_summary.txt"


def normalize_spacer(s):
    return str(s).strip().upper()


def safe_rank(series, ascending=False):
    return series.rank(method="min", ascending=ascending)


def main():
    print("Reading input files...")
    tool = pd.read_csv(INPUT_TOOL_TABLE)
    fotr = pd.read_csv(INPUT_FOTR_V2)

    # Normalize join keys
    tool["spacer_norm"] = tool["guide_sequence"].apply(normalize_spacer)
    fotr["spacer_norm"] = fotr["spacer"].apply(normalize_spacer)

    if "gene" in tool.columns:
        tool["gene_norm"] = tool["gene"].astype(str).str.strip()
    if "gene" in fotr.columns:
        fotr["gene_norm"] = fotr["gene"].astype(str).str.strip()

    # Select FOTR v2 columns defensively
    candidate_cols = [
        "gene_norm",
        "spacer_norm",
        "pam",
        "rna_accessibility_class",
        "spacer_accessibility_score",
        "target_region_bin",
        "functional_severity_v1",
        "target_context_weight_v1",
        "fotr_v1_priority_score",
        "fotr_v2_priority_score",
        "fotr_v2_rank",
        "fotr_v2_recommendation_status",
        "binding_risk",
        "rna_structure_risk",
        "risk_penalty",
    ]

    available_fotr_cols = [c for c in candidate_cols if c in fotr.columns]
    fotr_small = fotr[available_fotr_cols].drop_duplicates()

    print("Merging tool comparison with FOTR v2...")
    merged = tool.merge(
        fotr_small,
        on=["gene_norm", "spacer_norm"],
        how="left",
        suffixes=("", "_fotr")
    )

    # Build FOTR v2 rank if not present
    if "fotr_v2_rank" not in merged.columns and "fotr_v2_priority_score" in merged.columns:
        merged["fotr_v2_rank"] = merged.groupby("gene")["fotr_v2_priority_score"].rank(
            method="min", ascending=False
        )

    # Tool-rank columns available
    rank_cols = [c for c in ["rank_RS3", "rank_CFD", "rank_final", "fotr_v2_rank"] if c in merged.columns]

    # Disagreement metrics
    if rank_cols:
        merged["tool_rank_min"] = merged[rank_cols].min(axis=1)
        merged["tool_rank_max"] = merged[rank_cols].max(axis=1)
        merged["tool_rank_range"] = merged["tool_rank_max"] - merged["tool_rank_min"]
        merged["tool_rank_std"] = merged[rank_cols].std(axis=1)
    else:
        merged["tool_rank_range"] = np.nan
        merged["tool_rank_std"] = np.nan

    # Score spread across normalized scores where available
    score_cols = [c for c in ["RS3_score", "CFD_score", "MIT_score", "crisot_score", "final_score", "fotr_v2_priority_score"] if c in merged.columns]

    # Put scores on comparable 0-1 scale where possible
    normalized_score_cols = []
    for c in score_cols:
        new_c = c + "_norm_for_disagreement"
        vals = pd.to_numeric(merged[c], errors="coerce")
        if vals.max(skipna=True) > 1.5:
            merged[new_c] = vals / 100.0
        else:
            merged[new_c] = vals
        normalized_score_cols.append(new_c)

    if normalized_score_cols:
        merged["tool_score_range"] = merged[normalized_score_cols].max(axis=1) - merged[normalized_score_cols].min(axis=1)
        merged["tool_score_std"] = merged[normalized_score_cols].std(axis=1)
    else:
        merged["tool_score_range"] = np.nan
        merged["tool_score_std"] = np.nan

    # Combined disagreement index: rank spread + score spread, both normalized
    rank_range = pd.to_numeric(merged["tool_rank_range"], errors="coerce")
    score_range = pd.to_numeric(merged["tool_score_range"], errors="coerce")

    rank_component = rank_range / rank_range.max(skipna=True) if rank_range.max(skipna=True) not in [0, np.nan] else 0
    score_component = score_range / score_range.max(skipna=True) if score_range.max(skipna=True) not in [0, np.nan] else 0

    merged["tool_disagreement_index"] = (0.6 * rank_component.fillna(0)) + (0.4 * score_component.fillna(0))

    # Interpretation label
    def label_disagreement(x):
        if pd.isna(x):
            return "Not_available"
        if x >= 0.66:
            return "High_tool_disagreement"
        if x >= 0.33:
            return "Moderate_tool_disagreement"
        return "Low_tool_disagreement"

    merged["tool_disagreement_class"] = merged["tool_disagreement_index"].apply(label_disagreement)

    # Order useful columns first
    preferred_cols = [
        "gene", "position", "guide_sequence", "pam",
        "RS3_score", "CFD_score", "MIT_score", "crisot_score", "final_score", "fotr_v2_priority_score",
        "rank_RS3", "rank_CFD", "rank_final", "fotr_v2_rank",
        "tool_rank_range", "tool_score_range", "tool_disagreement_index", "tool_disagreement_class",
        "rna_accessibility_class", "spacer_accessibility_score",
        "target_region_bin", "functional_severity_v1", "target_context_weight_v1",
        "fotr_v2_recommendation_status"
    ]

    ordered = [c for c in preferred_cols if c in merged.columns]
    remaining = [c for c in merged.columns if c not in ordered and not c.endswith("_norm_for_disagreement")]
    final = merged[ordered + remaining]

    final.to_csv(OUT_TABLE, index=False)

    summary_lines = []
    summary_lines.append("Full Tool Comparison Summary")
    summary_lines.append("============================")
    summary_lines.append("")
    summary_lines.append(f"Input tool table: {INPUT_TOOL_TABLE}")
    summary_lines.append(f"Input FOTR v2 table: {INPUT_FOTR_V2}")
    summary_lines.append(f"Output table: {OUT_TABLE}")
    summary_lines.append("")
    summary_lines.append(f"Rows in final comparison table: {len(final)}")
    summary_lines.append(f"Unique spacers: {final['guide_sequence'].nunique() if 'guide_sequence' in final.columns else 'NA'}")
    summary_lines.append(f"Genes: {', '.join(sorted(final['gene'].dropna().astype(str).unique())) if 'gene' in final.columns else 'NA'}")
    summary_lines.append("")
    summary_lines.append("Disagreement class counts:")
    if "tool_disagreement_class" in final.columns:
        for k, v in final["tool_disagreement_class"].value_counts(dropna=False).items():
            summary_lines.append(f"- {k}: {v}")
    summary_lines.append("")
    summary_lines.append("Top 10 highest-disagreement guides:")
    top = final.sort_values("tool_disagreement_index", ascending=False).head(10)
    for _, row in top.iterrows():
        summary_lines.append(
            f"- {row.get('gene', 'NA')} | {row.get('guide_sequence', 'NA')} | "
            f"index={row.get('tool_disagreement_index', np.nan):.3f} | "
            f"class={row.get('tool_disagreement_class', 'NA')} | "
            f"FOTR_v2={row.get('fotr_v2_priority_score', 'NA')} | "
            f"RNA={row.get('rna_accessibility_class', 'NA')} | "
            f"recommendation={row.get('fotr_v2_recommendation_status', 'NA')}"
        )

    OUT_SUMMARY.write_text("\n".join(summary_lines))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(summary_lines[:20]))


if __name__ == "__main__":
    main()
