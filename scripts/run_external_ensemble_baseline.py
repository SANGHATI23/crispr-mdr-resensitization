#!/usr/bin/env python3
"""
External-style consensus / ensemble baseline for CRISPR-MDR resensitization.

Purpose:
Completes C. Consensus / ensemble baseline by creating a GuidePro-inspired
external ensemble-style comparator.

This does NOT claim direct GuidePro execution. It is an external-style
multi-source ensemble comparator for bacterial AMR guide prioritization.

Outputs:
  results_external/ensemble_baseline_comparison.csv
  results_external/top_guides_external_ensemble.csv
  results_external/ensemble_overlap_summary.csv
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


GUIDE_COL_CANDIDATES = [
    "guide_id", "sgRNA_id", "sgrna_id", "guide", "sgRNA", "sgrna",
    "spacer", "protospacer", "sequence", "guide_sequence"
]

INTERNAL_SCORE_CANDIDATES = [
    "final_score", "composite_score", "internal_consensus_score",
    "consensus_score", "priority_score", "internal_final_score"
]

RS3_SCORE_CANDIDATES = [
    "rs3_score", "Rule_Set_3_score", "rule_set_3_score",
    "ruleset3_score", "external_rs3_score"
]

SPECIFICITY_SCORE_CANDIDATES = [
    "specificity_score", "offtarget_specificity_score",
    "off_target_specificity_score", "cfd_specificity_score",
    "mit_specificity_score", "mit_score", "CFD_score", "cfd_score"
]

CONSERVATION_SCORE_CANDIDATES = [
    "conservation", "conservation_score", "strain_conservation_score",
    "conservation_proxy_score", "strain_aware_score", "strain_score"
]

TARGET_COL_CANDIDATES = [
    "gene", "target_gene", "amr_gene", "gene_name", "locus", "target"
]


def find_col(df: pd.DataFrame, candidates):
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def minmax(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(np.nan, index=series.index)

    min_v = x.min()
    max_v = x.max()

    if pd.isna(min_v) or pd.isna(max_v) or max_v == min_v:
        return pd.Series(0.5, index=series.index)

    z = (x - min_v) / (max_v - min_v)
    if not higher_is_better:
        z = 1.0 - z
    return z


def infer_specificity_direction(col_name: str) -> bool:
    """
    Returns True if higher score is better.
    For MIT score and specificity, higher is usually better.
    For raw off-target counts/risk columns, lower would be better.
    """
    name = col_name.lower()
    lower_is_better_terms = [
        "offtarget_count", "off_target_count", "mismatch_count",
        "risk", "penalty", "offtarget_risk", "off_target_risk"
    ]
    return not any(t in name for t in lower_is_better_terms)


def build_external_ensemble(df: pd.DataFrame) -> pd.DataFrame:
    guide_col = find_col(df, GUIDE_COL_CANDIDATES)
    internal_col = find_col(df, INTERNAL_SCORE_CANDIDATES)
    rs3_col = find_col(df, RS3_SCORE_CANDIDATES)
    specificity_col = find_col(df, SPECIFICITY_SCORE_CANDIDATES)
    conservation_col = find_col(df, CONSERVATION_SCORE_CANDIDATES)
    target_col = find_col(df, TARGET_COL_CANDIDATES)

    required = {
        "internal score": internal_col,
        "specificity score": specificity_col,
        "conservation score": conservation_col,
    }

    missing = [k for k, v in required.items() if v is None]
    if missing:
        print("\nERROR: Missing required score columns:", ", ".join(missing))
        print("\nAvailable columns:")
        for c in df.columns:
            print("  -", c)
        print("\nFix: rename your columns or add aliases inside this script.")
        sys.exit(1)

    out = df.copy()

    if guide_col is None:
        guide_col = "guide_id_auto"
        out[guide_col] = [f"guide_{i+1}" for i in range(len(out))]

    if target_col is None:
        target_col = "target_auto"
        out[target_col] = "unknown_target"

    out["norm_internal_score"] = minmax(out[internal_col], higher_is_better=True)
    out["norm_specificity_score"] = minmax(
        out[specificity_col],
        higher_is_better=infer_specificity_direction(specificity_col),
    )
    out["norm_conservation_score"] = minmax(out[conservation_col], higher_is_better=True)

    if rs3_col is not None and pd.to_numeric(out[rs3_col], errors="coerce").notna().sum() > 0:
        out["norm_rs3_score"] = minmax(out[rs3_col], higher_is_better=True)

        out["external_ensemble_score"] = (
            0.30 * out["norm_internal_score"]
            + 0.25 * out["norm_rs3_score"].fillna(out["norm_internal_score"])
            + 0.25 * out["norm_specificity_score"]
            + 0.20 * out["norm_conservation_score"]
        )
        out["ensemble_formula_used"] = (
            "0.30 internal + 0.25 RS3-style + 0.25 specificity + 0.20 conservation"
        )
        out["rs3_available"] = True
    else:
        out["norm_rs3_score"] = np.nan

        out["external_ensemble_score"] = (
            0.40 * out["norm_internal_score"]
            + 0.30 * out["norm_specificity_score"]
            + 0.30 * out["norm_conservation_score"]
        )
        out["ensemble_formula_used"] = (
            "0.40 internal + 0.30 specificity + 0.30 conservation"
        )
        out["rs3_available"] = False

    out["internal_consensus_rank"] = out["norm_internal_score"].rank(
        ascending=False, method="min"
    ).astype(int)

    out["external_ensemble_rank"] = out["external_ensemble_score"].rank(
        ascending=False, method="min"
    ).astype(int)

    out["rank_shift_external_minus_internal"] = (
        out["external_ensemble_rank"] - out["internal_consensus_rank"]
    )

    preferred_cols = [
        guide_col,
        target_col,
        internal_col,
        rs3_col,
        specificity_col,
        conservation_col,
        "norm_internal_score",
        "norm_rs3_score",
        "norm_specificity_score",
        "norm_conservation_score",
        "external_ensemble_score",
        "internal_consensus_rank",
        "external_ensemble_rank",
        "rank_shift_external_minus_internal",
        "rs3_available",
        "ensemble_formula_used",
    ]

    preferred_cols = [c for c in preferred_cols if c is not None and c in out.columns]
    remaining_cols = [c for c in out.columns if c not in preferred_cols]

    return out[preferred_cols + remaining_cols]


def overlap_summary(df: pd.DataFrame, ks=(10, 20, 50)) -> pd.DataFrame:
    rows = []

    for k in ks:
        internal_top = set(
            df.sort_values("internal_consensus_rank", ascending=True)
              .head(k)
              .index
        )
        external_top = set(
            df.sort_values("external_ensemble_rank", ascending=True)
              .head(k)
              .index
        )

        overlap = internal_top.intersection(external_top)
        jaccard = len(overlap) / len(internal_top.union(external_top)) if internal_top.union(external_top) else 0

        rows.append({
            "top_k": k,
            "internal_top_k_count": len(internal_top),
            "external_top_k_count": len(external_top),
            "overlap_count": len(overlap),
            "overlap_percent_of_k": round(100 * len(overlap) / k, 2),
            "jaccard_similarity": round(jaccard, 4),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Run external-style ensemble baseline for CRISPR-MDR guide ranking."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input master cross-model CSV file."
    )

    parser.add_argument(
        "--outdir",
        default="results_external",
        help="Output directory. Default: results_external"
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top guides to export. Default: 50"
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        sys.exit(1)

    df = pd.read_csv(input_path)

    if df.empty:
        print("ERROR: Input CSV is empty.")
        sys.exit(1)

    result = build_external_ensemble(df)

    comparison_path = outdir / "ensemble_baseline_comparison.csv"
    top_path = outdir / "top_guides_external_ensemble.csv"
    summary_path = outdir / "ensemble_overlap_summary.csv"

    result_sorted = result.sort_values(
        ["external_ensemble_rank", "internal_consensus_rank"],
        ascending=[True, True]
    )

    result_sorted.to_csv(comparison_path, index=False)
    result_sorted.head(args.top_n).to_csv(top_path, index=False)

    summary = overlap_summary(result)
    summary.to_csv(summary_path, index=False)

    print("\nExternal-style ensemble baseline complete.")
    print(f"Input rows: {len(df)}")
    print(f"Saved: {comparison_path}")
    print(f"Saved: {top_path}")
    print(f"Saved: {summary_path}")

    print("\nOverlap summary:")
    print(summary.to_string(index=False))

    print("\nSuggested paper wording:")
    print(
        "We implemented an external ensemble-style comparator inspired by "
        "published multi-source sgRNA prioritization frameworks. Because direct "
        "ensemble tools are mainly optimized for eukaryotic knockout contexts, "
        "the comparator was used as a ranking-control baseline rather than as "
        "a direct biological oracle for bacterial AMR targeting."
    )


if __name__ == "__main__":
    main()
