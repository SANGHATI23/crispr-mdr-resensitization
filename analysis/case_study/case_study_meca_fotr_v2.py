"""
mecA case study for FOTR v2.

Purpose:
- Create manuscript-ready mecA case-study outputs.
- Show why mecA is a strong biological example:
    1. Most mecA guides are structurally risky.
    2. The final top guide remains accessible.
    3. FOTR v2 keeps an interpretable candidate: accessible + early coding + recommended.
    4. Recommendation labels prevent structurally risky guides from being over-interpreted.

Input:
    results_fotr_v2/all_guides_fotr_v2_functional_context_recommended.csv

Outputs:
    results_case_study/meca_case_study_all_guides.csv
    results_case_study/meca_case_study_top_candidates.csv
    results_case_study/meca_case_study_summary.txt
"""

from pathlib import Path
import pandas as pd


INPUT_PATH = Path("results_fotr_v2/all_guides_fotr_v2_functional_context_recommended.csv")
OUTPUT_DIR = Path("results_case_study")

OUTPUT_ALL = OUTPUT_DIR / "meca_case_study_all_guides.csv"
OUTPUT_TOP = OUTPUT_DIR / "meca_case_study_top_candidates.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "meca_case_study_summary.txt"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    meca = df[df["gene"] == "mecA"].copy()

    if meca.empty:
        raise ValueError("No mecA rows found in input file.")

    meca = meca.sort_values("fotr_v2_priority_score", ascending=False)

    keep_cols = [
        "gene",
        "position",
        "strand",
        "spacer",
        "pam",
        "gc_content",
        "on_target_score",
        "specificity_score",
        "conservation_score",
        "final_score",
        "rna_accessibility_class",
        "seed_paired_fraction",
        "spacer_paired_fraction",
        "spacer_accessibility_score",
        "target_region_bin",
        "functional_severity_v1",
        "target_context_weight_v1",
        "target_functional_context_component_v2",
        "fotr_v1_priority_score",
        "fotr_v2_priority_score",
        "fotr_v1_rank_within_gene",
        "fotr_v2_rank_within_gene",
        "rank_change_v1_minus_v2",
        "fotr_v2_recommendation_status",
        "recommended_for_final_table",
        "crisot_available",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_mutant_score",
        "crisot_mean_delta_spec",
    ]

    existing_cols = [c for c in keep_cols if c in meca.columns]

    meca[existing_cols].to_csv(OUTPUT_ALL, index=False)

    top_candidates = (
        meca[meca["recommended_for_final_table"] == True]
        .sort_values("fotr_v2_priority_score", ascending=False)
        .head(5)
    )

    top_candidates[existing_cols].to_csv(OUTPUT_TOP, index=False)

    total_guides = len(meca)
    accessible_count = int((meca["rna_accessibility_class"] == "Accessible").sum())
    moderate_count = int((meca["rna_accessibility_class"] == "Moderate").sum())
    risky_count = int((meca["rna_accessibility_class"] == "Structurally_Risky").sum())

    accessible_pct = 100 * accessible_count / total_guides
    moderate_pct = 100 * moderate_count / total_guides
    risky_pct = 100 * risky_count / total_guides

    recommended_count = int(meca["recommended_for_final_table"].sum())
    high_priority_count = int(
        (meca["fotr_v2_recommendation_status"] == "High_priority_recommended").sum()
    )
    structurally_filtered_count = int(
        (meca["fotr_v2_recommendation_status"] == "Not_recommended_structural_risk").sum()
    )

    top = meca.iloc[0]

    top_original_rank = int(top["fotr_v1_rank_within_gene"])
    top_v2_rank = int(top["fotr_v2_rank_within_gene"])

    summary = f"""
mecA FOTR v2 Case Study Summary
================================

Total mecA guide rows: {total_guides}

RNA accessibility distribution:
- Accessible: {accessible_count} / {total_guides} ({accessible_pct:.2f}%)
- Moderate: {moderate_count} / {total_guides} ({moderate_pct:.2f}%)
- Structurally risky: {risky_count} / {total_guides} ({risky_pct:.2f}%)

Recommendation distribution:
- Final recommended candidates: {recommended_count}
- High-priority recommended candidates: {high_priority_count}
- Not recommended due to structural risk: {structurally_filtered_count}

Top mecA FOTR v2 candidate:
- Spacer: {top["spacer"]}
- PAM: {top["pam"]}
- Position: {top["position"]}
- RNA accessibility class: {top["rna_accessibility_class"]}
- Target region: {top["target_region_bin"]}
- Functional severity v1: {top["functional_severity_v1"]}
- Context weight v1: {top["target_context_weight_v1"]}
- FOTR v1 rank within mecA: {top_original_rank}
- FOTR v2 rank within mecA: {top_v2_rank}
- FOTR v2 priority score: {top["fotr_v2_priority_score"]}
- Recommendation status: {top["fotr_v2_recommendation_status"]}
- CRISOT available: {top["crisot_available"]}
- CRISOT WT specificity: {top["crisot_wt_specificity"]}

Manuscript-ready interpretation:
mecA provides a useful case study because RNA structural screening sharply narrows the candidate space. Among {total_guides} mecA guide rows, {risky_count} ({risky_pct:.2f}%) were classified as structurally risky, while only {accessible_count} ({accessible_pct:.2f}%) was classified as accessible. The final top FOTR v2 candidate, {top["spacer"]} with PAM {top["pam"]}, remained accessible, targeted an early coding region, and was labelled as {top["fotr_v2_recommendation_status"]}. This supports the use of RNA accessibility and functional-context filtering as interpretation layers beyond activity/conservation ranking alone.

Boundary statement:
This case study is computational prioritization only. RNAfold-derived accessibility is treated as a structural screening feature, not experimental proof of Cas9 loading, cleavage, antimicrobial resensitization, or clinical efficacy.
""".strip()

    OUTPUT_SUMMARY.write_text(summary)

    print(f"Wrote: {OUTPUT_ALL}")
    print(f"Wrote: {OUTPUT_TOP}")
    print(f"Wrote: {OUTPUT_SUMMARY}")

    print("\nmecA case-study summary:")
    print(summary)

    print("\nTop mecA recommended candidates:")
    print(top_candidates[existing_cols].to_string(index=False))


if __name__ == "__main__":
    main()