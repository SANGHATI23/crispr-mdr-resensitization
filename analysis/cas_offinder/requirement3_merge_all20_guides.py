#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

REQ2_ALL20 = PROJECT_ROOT / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2_database_refined/requirement2_database_refined_guide_summary_ALL20.csv"

REQ3_GUIDE = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/requirement3_guide_level_mobile_context_summary.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL20 = OUT_DIR / "requirement3_mobile_context_guide_summary_ALL20.csv"
OUT_TOP = OUT_DIR / "requirement3_top_guides_by_gene_ALL20.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_completion_summary.txt"


def pick_col(df, candidates, required=True):
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    if required:
        raise ValueError(
            f"Missing required column. Tried {candidates}. Available columns: {list(df.columns)}"
        )
    return None


def classify_final(row):
    total_hits = row.get("total_hits_req3", 0)
    mobile_burden = row.get("requirement3_mobile_context_burden", 0)
    req2_decision = str(row.get("decision", row.get("requirement2_decision", ""))).lower()

    if total_hits == 0 and mobile_burden == 0:
        return "Prioritize_zero_hit_low_mobile_context"

    if mobile_burden == 0:
        return "Prioritize_low_mobile_context"

    if mobile_burden <= 10:
        return "Backup_moderate_mobile_context"

    if mobile_burden <= 50:
        return "Caution_elevated_mobile_context"

    return "Deprioritize_high_mobile_context"


def main():
    print("Requirement 3 Step 2: merging mobile-context burden into all-20 guide table")
    print(f"Project root: {PROJECT_ROOT}")

    if not REQ2_ALL20.exists():
        raise FileNotFoundError(f"Requirement 2 all-20 guide table not found: {REQ2_ALL20}")

    if not REQ3_GUIDE.exists():
        raise FileNotFoundError(f"Requirement 3 guide summary not found: {REQ3_GUIDE}")

    req2 = pd.read_csv(REQ2_ALL20)
    req3 = pd.read_csv(REQ3_GUIDE)

    print(f"Loaded Requirement 2 all-20 table: {req2.shape}")
    print(f"Loaded Requirement 3 guide table: {req3.shape}")

    req2_gene_col = pick_col(req2, ["target_gene", "gene", "amr_gene"])
    req2_guide_col = pick_col(req2, ["guide_id", "selected_guide_id", "guide"])

    req3_gene_col = pick_col(req3, ["target_gene", "gene", "amr_gene"])
    req3_guide_col = pick_col(req3, ["guide_id", "selected_guide_id", "guide"])

    # Standardize merge column names.
    req2 = req2.rename(columns={req2_gene_col: "target_gene", req2_guide_col: "guide_id"})
    req3 = req3.rename(columns={req3_gene_col: "target_gene", req3_guide_col: "guide_id"})

    if "total_hits" in req3.columns:
        req3 = req3.rename(columns={"total_hits": "total_hits_req3"})

    merged = req2.merge(
        req3,
        on=["target_gene", "guide_id"],
        how="left",
        suffixes=("_req2", "_req3")
    )

    # Fill zero-hit guide after merge.
    fill_zero_cols = [
        "total_hits_req3",
        "mobile_context_hits",
        "plasmid_context_hits",
        "integron_context_hits",
        "sccmec_context_hits",
        "resistance_island_like_hits",
        "operon_proxy_hits",
        "requirement3_mobile_context_burden",
    ]

    for col in fill_zero_cols:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    if "requirement3_mobile_context_decision" in merged.columns:
        merged["requirement3_mobile_context_decision"] = merged[
            "requirement3_mobile_context_decision"
        ].fillna("Low_mobile_context_signal_zero_hit")

    merged["requirement3_final_decision"] = merged.apply(classify_final, axis=1)

    # Ranking logic:
    # lower mobile burden is better;
    # fewer Req3 hits is better;
    # then keep any previous burden/rank if available.
    sort_cols = ["target_gene", "requirement3_mobile_context_burden", "total_hits_req3"]

    for optional in [
        "db_refined_burden",
        "database_refined_burden",
        "functional_offtarget_burden",
        "total_burden",
        "burden",
    ]:
        if optional in merged.columns:
            sort_cols.append(optional)
            break

    merged = merged.sort_values(sort_cols, ascending=True).copy()

    merged["requirement3_rank_within_gene"] = (
        merged.groupby("target_gene").cumcount() + 1
    )

    top_guides = (
        merged.sort_values(
            ["target_gene", "requirement3_rank_within_gene"],
            ascending=True
        )
        .groupby("target_gene")
        .head(2)
        .copy()
    )

    merged.to_csv(OUT_ALL20, index=False)
    top_guides.to_csv(OUT_TOP, index=False)

    expected_guides = 20
    observed_guides = merged["guide_id"].nunique()
    zero_hit_guides = int((merged["total_hits_req3"] == 0).sum())

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Completion Summary\n")
        f.write("================================\n\n")
        f.write("Step completed: Requirement 3 mobile-context guide-level all-20 merge\n\n")
        f.write(f"Requirement 2 all-20 input: {REQ2_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Requirement 3 guide summary input: {REQ3_GUIDE.relative_to(PROJECT_ROOT)}\n\n")
        f.write(f"Expected selected guides: {expected_guides}\n")
        f.write(f"Observed unique guides after merge: {observed_guides}\n")
        f.write(f"Zero-hit guides retained: {zero_hit_guides}\n\n")

        f.write("Requirement 3 decision counts:\n")
        for k, v in merged["requirement3_final_decision"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nTop 2 guides per gene after Requirement 3 mobile-context mapping:\n")
        keep_cols = [
            "target_gene",
            "requirement3_rank_within_gene",
            "guide_id",
            "total_hits_req3",
            "mobile_context_hits",
            "plasmid_context_hits",
            "integron_context_hits",
            "sccmec_context_hits",
            "resistance_island_like_hits",
            "requirement3_mobile_context_burden",
            "requirement3_final_decision",
        ]
        keep_cols = [c for c in keep_cols if c in top_guides.columns]

        for _, row in top_guides[keep_cols].iterrows():
            f.write(
                f"- {row['target_gene']} rank {row['requirement3_rank_within_gene']}: "
                f"{row['guide_id']} | hits={row['total_hits_req3']} | "
                f"mobile_burden={row['requirement3_mobile_context_burden']} | "
                f"decision={row['requirement3_final_decision']}\n"
            )

        f.write("\nOutput files:\n")
        f.write(f"- {OUT_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_TOP.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_SUMMARY.relative_to(PROJECT_ROOT)}\n")

    print(f"Wrote: {OUT_ALL20.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_TOP.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nGuide count check:")
    print(f"Expected guides: {expected_guides}")
    print(f"Observed unique guides: {observed_guides}")
    print(f"Zero-hit guides retained: {zero_hit_guides}")

    print("\nTop 2 guides per gene after Requirement 3:")
    cols = [
        "target_gene",
        "requirement3_rank_within_gene",
        "guide_id",
        "total_hits_req3",
        "mobile_context_hits",
        "plasmid_context_hits",
        "integron_context_hits",
        "sccmec_context_hits",
        "resistance_island_like_hits",
        "requirement3_mobile_context_burden",
        "requirement3_final_decision",
    ]
    cols = [c for c in cols if c in top_guides.columns]
    print(top_guides[cols].to_string(index=False))

    print("\nRequirement 3 Step 2 completed.")


if __name__ == "__main__":
    main()
