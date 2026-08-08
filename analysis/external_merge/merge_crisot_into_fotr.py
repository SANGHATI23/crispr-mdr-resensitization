from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

FOTR_FILE = PROJECT_ROOT / "results_fotr" / "all_guides_fotr_scores.csv"

CRISOT_FILE = (
    PROJECT_ROOT
    / "results_external"
    / "crisot"
    / "merged"
    / "crisot_guide_level_summary.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "results_fotr"
    / "all_guides_fotr_scores_with_crisot.csv"
)


def main() -> None:
    fotr = pd.read_csv(FOTR_FILE)
    crisot = pd.read_csv(CRISOT_FILE)

    fotr["spacer"] = fotr["spacer"].astype(str).str.upper().str.strip()
    fotr["pam"] = fotr["pam"].astype(str).str.upper().str.strip()
    fotr["guide_23mer"] = fotr["spacer"] + fotr["pam"]

    crisot["spacer"] = crisot["spacer"].astype(str).str.upper().str.strip()
    crisot["pam"] = crisot["pam"].astype(str).str.upper().str.strip()
    crisot["crisot_guide_23mer"] = crisot["guide_23mer"].astype(str).str.upper().str.strip()
    crisot["crisot_pam_used"] = crisot["pam"]

    crisot_cols = [
        "spacer",
        "crisot_guide_23mer",
        "crisot_pam_used",
        "crisot_wt_score",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_mutant_score",
        "crisot_min_mutant_score",
        "crisot_max_mutant_score",
        "crisot_mean_delta_spec",
        "crisot_max_delta_spec",
    ]

    merged = fotr.merge(
        crisot[crisot_cols],
        on="spacer",
        how="left",
    )

    merged["crisot_available"] = merged["crisot_wt_score"].notna()
    merged["crisot_merge_level"] = merged["crisot_available"].map(
        {True: "spacer_level", False: "not_available"}
    )

    merged.to_csv(OUTPUT_FILE, index=False)

    print(f"FOTR input guides: {len(fotr)}")
    print(f"CRISOT guide summaries: {len(crisot)}")
    print(f"Merged output guides: {len(merged)}")
    print(f"CRISOT available: {merged['crisot_available'].sum()}")
    print(f"CRISOT missing: {(~merged['crisot_available']).sum()}")
    print(f"Output saved: {OUTPUT_FILE}")

    print()
    print("Top rows with CRISOT spacer-level merge:")
    cols = [
        "gene",
        "spacer",
        "pam",
        "fotr_priority_score",
        "rna_accessibility_class",
        "crisot_available",
        "crisot_merge_level",
        "crisot_pam_used",
        "crisot_wt_score",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_delta_spec",
    ]

    print(merged[cols].head(12).to_string(index=False))


if __name__ == "__main__":
    main()