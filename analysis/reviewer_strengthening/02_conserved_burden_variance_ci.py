
import pandas as pd
import numpy as np
from pathlib import Path

OUT_DIR = Path("results_reviewer_strengthening/variance_ci")
OUT_DIR.mkdir(parents=True, exist_ok=True)

HIT_FILE = Path("results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/requirement4_hit_level_conservation_ready_table.csv")
GUIDE_FILE = Path("results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/final_requirement4_complete/requirement4_final_integrated_conserved_offtarget_table_ALL20.csv")

def main():
    hits = pd.read_csv(HIT_FILE)
    guides = pd.read_csv(GUIDE_FILE)

    # Ensure numeric burden
    hits["db_refined_burden_score"] = pd.to_numeric(
        hits["db_refined_burden_score"], errors="coerce"
    ).fillna(0)

    # Total accessions screened from guide table
    total_accessions = int(guides["total_accessions_screened"].max())

    # Build guide-accession burden table for accessions WITH hits
    burden_hit_accessions = (
        hits.groupby(["target_gene", "guide_id", "spacer", "accession"])
        .agg(
            accession_hit_count=("off_target_sequence", "count"),
            accession_burden_sum=("db_refined_burden_score", "sum"),
            accession_burden_mean=("db_refined_burden_score", "mean"),
            accession_high_confidence_hits=("db_refined_high_confidence_functional_hit", lambda x: x.astype(str).str.lower().isin(["true", "1", "yes"]).sum())
        )
        .reset_index()
    )

    # Include zero-burden accessions for guides with no hit in some genomes.
    # We do not know all accession names for zero-hit guide/accession pairs from hit table,
    # so we use observed accession list from the full hit-level table as the screened accession set.
    accession_list = sorted(hits["accession"].dropna().unique().tolist())

    guide_list = guides[["target_gene", "guide_id", "spacer"]].drop_duplicates().copy()
    grid = guide_list.assign(key=1).merge(
        pd.DataFrame({"accession": accession_list, "key": 1}),
        on="key"
    ).drop(columns=["key"])

    burden_all_accessions = grid.merge(
        burden_hit_accessions,
        on=["target_gene", "guide_id", "spacer", "accession"],
        how="left"
    )

    for col in ["accession_hit_count", "accession_burden_sum", "accession_burden_mean", "accession_high_confidence_hits"]:
        burden_all_accessions[col] = burden_all_accessions[col].fillna(0)

    # Per-guide variance/CI
    rows = []
    for (target_gene, guide_id, spacer), g in burden_all_accessions.groupby(["target_gene", "guide_id", "spacer"]):
        n = len(g)
        mean_burden = g["accession_burden_sum"].mean()
        sd_burden = g["accession_burden_sum"].std(ddof=1) if n > 1 else 0
        se_burden = sd_burden / np.sqrt(n) if n > 0 else 0
        ci95_low = mean_burden - 1.96 * se_burden
        ci95_high = mean_burden + 1.96 * se_burden

        rows.append({
            "target_gene": target_gene,
            "guide_id": guide_id,
            "spacer": spacer,
            "n_accessions_used_for_ci": n,
            "mean_accession_burden": round(mean_burden, 4),
            "sd_accession_burden": round(sd_burden, 4),
            "se_accession_burden": round(se_burden, 4),
            "ci95_low": round(max(0, ci95_low), 4),
            "ci95_high": round(ci95_high, 4),
            "accessions_with_nonzero_burden": int((g["accession_burden_sum"] > 0).sum()),
            "fraction_accessions_with_nonzero_burden": round((g["accession_burden_sum"] > 0).mean(), 4),
            "max_accession_burden": round(g["accession_burden_sum"].max(), 4),
            "total_accession_burden_sum": round(g["accession_burden_sum"].sum(), 4),
        })

    ci = pd.DataFrame(rows)

    # Add original integrated conserved burden columns
    merge_cols = [
        "target_gene",
        "guide_id",
        "conserved_functional_burden_score",
        "total_offtarget_hit_rows",
        "accessions_with_any_hit",
        "fraction_accessions_with_any_hit",
        "final_req4_integrated_decision",
        "final_req4_rank_within_gene"
    ]

    merge_cols_existing = [c for c in merge_cols if c in guides.columns]
    ci = ci.merge(guides[merge_cols_existing], on=["target_gene", "guide_id"], how="left")

    ci = ci.sort_values(["target_gene", "final_req4_rank_within_gene", "mean_accession_burden"], na_position="last")

    ci_path = OUT_DIR / "conserved_burden_variance_ci_by_guide.csv"
    acc_path = OUT_DIR / "guide_accession_burden_matrix.csv"
    ci.to_csv(ci_path, index=False)
    burden_all_accessions.to_csv(acc_path, index=False)

    # Gene-level summary
    gene_summary = (
        ci.groupby("target_gene")
        .agg(
            guides=("guide_id", "nunique"),
            mean_of_mean_accession_burden=("mean_accession_burden", "mean"),
            max_mean_accession_burden=("mean_accession_burden", "max"),
            mean_fraction_accessions_nonzero=("fraction_accessions_with_nonzero_burden", "mean"),
            max_fraction_accessions_nonzero=("fraction_accessions_with_nonzero_burden", "max"),
        )
        .reset_index()
    )

    gene_path = OUT_DIR / "conserved_burden_variance_ci_gene_summary.csv"
    gene_summary.to_csv(gene_path, index=False)

    # Text summary
    txt = []
    txt.append("Conserved off-target burden variance and confidence interval analysis")
    txt.append("====================================================================")
    txt.append("")
    txt.append(f"Input hit-level table: {HIT_FILE}")
    txt.append(f"Input final guide-level table: {GUIDE_FILE}")
    txt.append(f"Guides evaluated: {ci['guide_id'].nunique()}")
    txt.append(f"Accessions used for CI grid: {len(accession_list)}")
    txt.append("")
    txt.append("Purpose:")
    txt.append("This analysis estimates guide-level variability of database-refined off-target burden across screened accessions.")
    txt.append("It reports mean accession burden, standard deviation, standard error, and approximate 95% confidence intervals.")
    txt.append("")
    txt.append("Reviewer-facing interpretation:")
    txt.append("This addresses whether final conserved-burden scores are driven by broad accession-level recurrence or by a few high-burden accessions.")
    txt.append("Lower mean burden and narrower intervals support cleaner computational candidates.")
    txt.append("")
    txt.append("Top guide-level CI table:")
    display_cols = [
        "target_gene", "guide_id", "mean_accession_burden", "sd_accession_burden",
        "ci95_low", "ci95_high", "fraction_accessions_with_nonzero_burden",
        "conserved_functional_burden_score", "final_req4_integrated_decision"
    ]
    display_cols_existing = [c for c in display_cols if c in ci.columns]
    txt.append(ci[display_cols_existing].to_string(index=False))
    txt.append("")
    txt.append("Recommended manuscript wording:")
    txt.append("To evaluate whether conserved off-target burden was concentrated in a few genomes or distributed across the accession panel, we computed per-guide accession-level burden statistics. For each guide, database-refined burden was aggregated within each screened accession and summarized using mean, standard deviation, standard error, and approximate 95% confidence intervals. This analysis provides a stability check for the final integrated ranking by distinguishing consistently low-burden candidates from guides whose burden is driven by recurrent or high-intensity accession-level signals.")

    txt_path = OUT_DIR / "conserved_burden_variance_ci_summary.txt"
    txt_path.write_text("\n".join(txt))

    print("Wrote:")
    print(ci_path)
    print(acc_path)
    print(gene_path)
    print(txt_path)
    print()
    print(ci[[
        "target_gene", "guide_id", "mean_accession_burden", "sd_accession_burden",
        "ci95_low", "ci95_high", "fraction_accessions_with_nonzero_burden"
    ]].to_string(index=False))

if __name__ == "__main__":
    main()
