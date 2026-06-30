#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3/requirement3_mobile_context_guide_summary_ALL20.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_combined"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL20 = OUT_DIR / "requirement3_combined_requirement2_requirement3_final_guide_table_ALL20.csv"
OUT_TOP = OUT_DIR / "requirement3_combined_top_guides_by_gene.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_combined_final_summary.txt"


def combined_decision(row):
    req2 = str(row.get("db_refined_decision", "")).upper()
    req3_burden = float(row.get("requirement3_mobile_context_burden", 0))
    total_hits = float(row.get("total_hits_req3", 0))

    if total_hits == 0 and req3_burden == 0:
        return "FINAL_PRIORITIZE_ZERO_HIT_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT"

    if "DEPRIORITIZE" in req2 and req3_burden > 50:
        return "FINAL_DEPRIORITIZE_HIGH_FUNCTIONAL_AND_HIGH_MOBILE_CONTEXT"

    if "DEPRIORITIZE" in req2:
        return "FINAL_DEPRIORITIZE_HIGH_FUNCTIONAL_BURDEN"

    if req3_burden > 50:
        return "FINAL_DEPRIORITIZE_HIGH_MOBILE_CONTEXT"

    if "USE_WITH_CAUTION" in req2 or req3_burden > 10:
        return "FINAL_CAUTION_FUNCTIONAL_OR_MOBILE_CONTEXT"

    if "KEEP_AS_BACKUP" in req2 or req3_burden > 0:
        return "FINAL_BACKUP_ACCEPTABLE_WITH_MONITORING"

    if "PRIORITIZE" in req2 and req3_burden == 0:
        return "FINAL_PRIORITIZE_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT"

    return "FINAL_BACKUP_REVIEW"


def decision_order(decision):
    order = {
        "FINAL_PRIORITIZE_ZERO_HIT_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT": 0,
        "FINAL_PRIORITIZE_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT": 1,
        "FINAL_BACKUP_ACCEPTABLE_WITH_MONITORING": 2,
        "FINAL_CAUTION_FUNCTIONAL_OR_MOBILE_CONTEXT": 3,
        "FINAL_BACKUP_REVIEW": 4,
        "FINAL_DEPRIORITIZE_HIGH_FUNCTIONAL_BURDEN": 5,
        "FINAL_DEPRIORITIZE_HIGH_MOBILE_CONTEXT": 6,
        "FINAL_DEPRIORITIZE_HIGH_FUNCTIONAL_AND_HIGH_MOBILE_CONTEXT": 7,
    }
    return order.get(decision, 9)


def main():
    print("Creating combined Requirement 2 + Requirement 3 final guide table...")

    if not INPUT.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT}")

    df = pd.read_csv(INPUT)
    print(f"Loaded input: {df.shape}")

    df["combined_req2_req3_final_decision"] = df.apply(combined_decision, axis=1)
    df["combined_req2_req3_decision_order"] = df["combined_req2_req3_final_decision"].apply(decision_order)

    sort_cols = [
        "target_gene",
        "combined_req2_req3_decision_order",
        "requirement3_mobile_context_burden",
    ]

    if "db_refined_total_burden" in df.columns:
        sort_cols.append("db_refined_total_burden")

    if "total_hits_req3" in df.columns:
        sort_cols.append("total_hits_req3")

    df = df.sort_values(sort_cols, ascending=True).copy()

    df["combined_req2_req3_rank_within_gene"] = df.groupby("target_gene").cumcount() + 1

    top = df.groupby("target_gene").head(2).copy()

    df.to_csv(OUT_ALL20, index=False)
    top.to_csv(OUT_TOP, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Combined Requirement 2 + Requirement 3 Final Summary\n")
        f.write("====================================================\n\n")

        f.write("Purpose:\n")
        f.write(
            "This table combines database-refined functional off-target burden "
            "from Requirement 2 with plasmid/mobile/integron/SCCmec/resistance-island "
            "context burden from Requirement 3.\n\n"
        )

        f.write(f"Input table: {INPUT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output all-20 table: {OUT_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output top-guide table: {OUT_TOP.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Total guide rows: {len(df)}\n")
        f.write(f"Unique guides: {df['guide_id'].nunique()}\n")
        f.write(f"Target genes: {df['target_gene'].nunique()}\n\n")

        f.write("Combined decision counts:\n")
        for k, v in df["combined_req2_req3_final_decision"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nFinal top 2 guides per gene after combined Requirement 2 + 3 scoring:\n")
        for _, r in top.iterrows():
            f.write(
                f"- {r['target_gene']} rank {r['combined_req2_req3_rank_within_gene']}: "
                f"{r['guide_id']} | "
                f"Req2_decision={r.get('db_refined_decision', '')} | "
                f"Req2_burden={r.get('db_refined_total_burden', '')} | "
                f"Req3_mobile_burden={r.get('requirement3_mobile_context_burden', '')} | "
                f"Final={r['combined_req2_req3_final_decision']}\n"
            )

    print(f"Wrote: {OUT_ALL20.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_TOP.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nFinal top 2 guides per gene:")
    show_cols = [
        "target_gene",
        "combined_req2_req3_rank_within_gene",
        "guide_id",
        "total_hits_req3",
        "db_refined_total_burden",
        "db_refined_decision",
        "requirement3_mobile_context_burden",
        "combined_req2_req3_final_decision",
    ]
    show_cols = [c for c in show_cols if c in top.columns]
    print(top[show_cols].to_string(index=False))

    print("\nCombined Requirement 2 + 3 final decision completed.")


if __name__ == "__main__":
    main()
