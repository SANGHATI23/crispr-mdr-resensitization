#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path("/Users/sanghati/Documents/crispr-mdr-resensitization")

REQ1_GUIDES = PROJECT_ROOT / "results_cas_offinder/expanded_panel/final_requirement1/requirement1_expanded_casoffinder_guide_level_summary.csv"

REQ2_GUIDES = PROJECT_ROOT / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2/requirement2_final_guide_prioritization_table.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2"

OUT_FIXED = OUT_DIR / "requirement2_final_guide_prioritization_table_ALL20.csv"
OUT_TOP = OUT_DIR / "requirement2_top_guides_by_gene_ALL20.csv"
OUT_ZERO = OUT_DIR / "requirement2_zero_hit_guides_retained.csv"
OUT_SUMMARY = OUT_DIR / "requirement2_ALL20_guide_audit_summary.txt"

req1 = pd.read_csv(REQ1_GUIDES)
req2 = pd.read_csv(REQ2_GUIDES)

print("Requirement 1 guide rows:", req1.shape)
print("Requirement 2 guide rows before fix:", req2.shape)

# Identify guide/gene columns
guide_col = "guide_id"
gene_col = "target_gene"

if guide_col not in req1.columns or guide_col not in req2.columns:
    raise ValueError("guide_id column missing from one of the files.")

if gene_col not in req1.columns or gene_col not in req2.columns:
    raise ValueError("target_gene column missing from one of the files.")

# Keep all Requirement 1 guide rows, merge Requirement 2 burden columns
merged = req1[[gene_col, guide_col]].drop_duplicates().merge(
    req2,
    on=[gene_col, guide_col],
    how="left",
    suffixes=("", "_req2")
)

numeric_fill_cols = [
    "total_offtarget_hits",
    "coding_cds_hits",
    "noncoding_or_intergenic_hits",
    "amr_hits",
    "virulence_hits",
    "essential_like_hits",
    "mobile_or_plasmid_hits",
    "low_mismatch_hits",
    "high_interest_hits",
    "total_functional_burden",
    "mean_functional_burden",
]

for col in numeric_fill_cols:
    if col in merged.columns:
        merged[col] = merged[col].fillna(0)

# Identify zero-hit guides that were absent from Requirement 2 groupby
missing_mask = merged["functional_annotation_recommendation"].isna()

merged.loc[missing_mask, "functional_annotation_recommendation"] = (
    "Strong_candidate_zero_detected_offtarget_burden"
)
merged.loc[missing_mask, "recommendation_clean"] = (
    "Strong zero-hit candidate"
)
merged.loc[missing_mask, "requirement2_decision"] = (
    "PRIORITIZE_OR_KEEP_AS_LOW_OFFTARGET_CANDIDATE"
)
merged.loc[missing_mask, "recommendation_order"] = 0

# Re-rank after adding zero-hit guide
merged = merged.sort_values(
    by=[
        gene_col,
        "recommendation_order",
        "total_functional_burden",
        "high_interest_hits",
        "total_offtarget_hits",
    ],
    ascending=True
).copy()

merged["requirement2_rank_within_gene"] = merged.groupby(gene_col).cumcount() + 1

top = merged[merged["requirement2_rank_within_gene"] <= 2].copy()
zero = merged[merged["functional_annotation_recommendation"] == "Strong_candidate_zero_detected_offtarget_burden"].copy()

merged.to_csv(OUT_FIXED, index=False)
top.to_csv(OUT_TOP, index=False)
zero.to_csv(OUT_ZERO, index=False)

with open(OUT_SUMMARY, "w") as f:
    f.write("Requirement 2 ALL-20 Guide Audit Summary\n")
    f.write("=======================================\n\n")
    f.write(f"Requirement 1 guide rows expected: {len(req1)}\n")
    f.write(f"Requirement 2 guide rows before zero-hit fix: {len(req2)}\n")
    f.write(f"Requirement 2 guide rows after zero-hit fix: {len(merged)}\n")
    f.write(f"Zero-hit guides retained: {len(zero)}\n\n")
    f.write("Zero-hit guides:\n")
    if len(zero) > 0:
        f.write(zero[[gene_col, guide_col, "requirement2_decision"]].to_string(index=False))
    else:
        f.write("None\n")
    f.write("\n\nTop guides by gene after retaining zero-hit guides:\n")
    f.write(top[[gene_col, guide_col, "requirement2_rank_within_gene", "total_offtarget_hits", "high_interest_hits", "total_functional_burden", "recommendation_clean", "requirement2_decision"]].to_string(index=False))
    f.write("\n\nQC status: PASS if after-fix guide rows = 20.\n")

print("\nWrote:")
print(OUT_FIXED)
print(OUT_TOP)
print(OUT_ZERO)
print(OUT_SUMMARY)

print("\nAfter-fix guide rows:", len(merged))
print("\nZero-hit retained:")
print(zero[[gene_col, guide_col, "requirement2_decision"]].to_string(index=False))

print("\nTop guides by gene after ALL20 fix:")
print(top[[gene_col, guide_col, "requirement2_rank_within_gene", "total_offtarget_hits", "high_interest_hits", "total_functional_burden", "recommendation_clean", "requirement2_decision"]].to_string(index=False))
