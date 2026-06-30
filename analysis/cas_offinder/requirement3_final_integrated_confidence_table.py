#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_integron_is_transposon_mapping/requirement3_all20_with_direct_plasmid_sccmec_integron_is_transposon_evidence.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_complete"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL20 = OUT_DIR / "requirement3_final_integrated_mobile_element_mapping_ALL20.csv"
OUT_TOP = OUT_DIR / "requirement3_final_integrated_top_guides_by_gene.csv"
OUT_DEPRIORITIZED = OUT_DIR / "requirement3_final_deprioritized_mobile_context_guides.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_final_integrated_completion_summary.txt"


def get_num(row, col):
    try:
        return float(row.get(col, 0))
    except Exception:
        return 0.0


def evidence_score(row):
    score = 0

    # Existing proxy/context burden.
    req3_burden = get_num(row, "requirement3_mobile_context_burden")

    if req3_burden > 50:
        score += 3
    elif req3_burden > 10:
        score += 2
    elif req3_burden > 0:
        score += 1

    # Direct plasmid evidence.
    plasmid_strength = str(row.get("direct_plasmid_evidence_strength", ""))
    if "Strong" in plasmid_strength:
        score += 3
    elif "Moderate" in plasmid_strength:
        score += 2
    elif "Low" in plasmid_strength:
        score += 1

    # Direct SCCmec evidence.
    sccmec_strength = str(row.get("direct_sccmec_evidence_strength", ""))
    if "Strong" in sccmec_strength:
        score += 3
    elif "Moderate" in sccmec_strength:
        score += 2
    elif "Low" in sccmec_strength:
        score += 1

    # Direct mobile element evidence.
    mobile_strength = str(row.get("direct_mobile_element_evidence_strength", ""))
    if "Strong" in mobile_strength:
        score += 3
    elif "Moderate" in mobile_strength:
        score += 2
    elif "Low" in mobile_strength:
        score += 1

    # Extra penalty for high functional burden from Requirement 2.
    req2_decision = str(row.get("db_refined_decision", "")).upper()
    if "DEPRIORITIZE" in req2_decision:
        score += 3
    elif "USE_WITH_CAUTION" in req2_decision:
        score += 2
    elif "KEEP_AS_BACKUP" in req2_decision:
        score += 1

    return score


def confidence_label(score):
    if score == 0:
        return "No_detected_mobile_context_or_functional_burden"
    if score <= 2:
        return "Low_integrated_mobile_context_confidence"
    if score <= 5:
        return "Moderate_integrated_mobile_context_confidence"
    if score <= 8:
        return "High_integrated_mobile_context_confidence"
    return "Very_high_integrated_mobile_context_confidence"


def final_decision(row):
    score = row["requirement3_integrated_evidence_score"]
    total_hits = get_num(row, "total_hits_req3")
    req2_decision = str(row.get("db_refined_decision", "")).upper()

    if total_hits == 0 and score == 0:
        return "FINAL_PRIORITIZE_ZERO_HIT_NO_MOBILE_CONTEXT_EVIDENCE"

    if score == 0 and "PRIORITIZE" in req2_decision:
        return "FINAL_PRIORITIZE_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT"

    if score <= 2 and "DEPRIORITIZE" not in req2_decision:
        return "FINAL_BACKUP_OR_PRIORITIZE_WITH_LOW_CONTEXT_RISK"

    if score <= 5:
        return "FINAL_CAUTION_MODERATE_CONTEXT_OR_FUNCTIONAL_RISK"

    if score <= 8:
        return "FINAL_DEPRIORITIZE_HIGH_CONTEXT_OR_FUNCTIONAL_RISK"

    return "FINAL_DEPRIORITIZE_VERY_HIGH_CONTEXT_AND_FUNCTIONAL_RISK"


def main():
    print("Creating final integrated Requirement 3 confidence table...")

    if not INPUT.exists():
        raise FileNotFoundError(f"Input table not found: {INPUT}")

    df = pd.read_csv(INPUT)
    print(f"Loaded input: {df.shape}")

    df["requirement3_integrated_evidence_score"] = df.apply(evidence_score, axis=1)
    df["requirement3_integrated_confidence_label"] = df[
        "requirement3_integrated_evidence_score"
    ].apply(confidence_label)

    df["requirement3_final_integrated_decision"] = df.apply(final_decision, axis=1)

    decision_order = {
        "FINAL_PRIORITIZE_ZERO_HIT_NO_MOBILE_CONTEXT_EVIDENCE": 0,
        "FINAL_PRIORITIZE_LOW_FUNCTIONAL_AND_MOBILE_CONTEXT": 1,
        "FINAL_BACKUP_OR_PRIORITIZE_WITH_LOW_CONTEXT_RISK": 2,
        "FINAL_CAUTION_MODERATE_CONTEXT_OR_FUNCTIONAL_RISK": 3,
        "FINAL_DEPRIORITIZE_HIGH_CONTEXT_OR_FUNCTIONAL_RISK": 4,
        "FINAL_DEPRIORITIZE_VERY_HIGH_CONTEXT_AND_FUNCTIONAL_RISK": 5,
    }

    df["requirement3_final_decision_order"] = df[
        "requirement3_final_integrated_decision"
    ].map(decision_order).fillna(9).astype(int)

    sort_cols = [
        "target_gene",
        "requirement3_final_decision_order",
        "requirement3_integrated_evidence_score",
    ]

    if "db_refined_total_burden" in df.columns:
        sort_cols.append("db_refined_total_burden")

    if "total_hits_req3" in df.columns:
        sort_cols.append("total_hits_req3")

    df = df.sort_values(sort_cols, ascending=True).copy()
    df["requirement3_final_integrated_rank_within_gene"] = (
        df.groupby("target_gene").cumcount() + 1
    )

    top = df.groupby("target_gene").head(2).copy()

    deprioritized = df[
        df["requirement3_final_integrated_decision"].str.contains("DEPRIORITIZE", na=False)
    ].copy()

    df.to_csv(OUT_ALL20, index=False)
    top.to_csv(OUT_TOP, index=False)
    deprioritized.to_csv(OUT_DEPRIORITIZED, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Final Integrated Completion Summary\n")
        f.write("================================================\n\n")

        f.write("Purpose:\n")
        f.write(
            "This final table integrates staged GFF-neighborhood mobile-context mapping, "
            "direct plasmid-contig evidence, direct SCCmec-like evidence, direct integron/"
            "insertion-sequence/transposon evidence, and Requirement 2 database-refined "
            "functional burden into one all-20 guide-level decision table.\n\n"
        )

        f.write(f"Input table: {INPUT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output all-20 final table: {OUT_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output top-guide table: {OUT_TOP.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Output deprioritized-guide table: {OUT_DEPRIORITIZED.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Final guide rows: {len(df)}\n")
        f.write(f"Unique guides: {df['guide_id'].nunique()}\n")
        f.write(f"Target genes represented: {df['target_gene'].nunique()}\n")
        f.write(f"Zero-hit guides retained: {int((df.get('total_hits_req3', pd.Series([0]*len(df))) == 0).sum())}\n\n")

        f.write("Integrated confidence label counts:\n")
        for k, v in df["requirement3_integrated_confidence_label"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nFinal integrated decision counts:\n")
        for k, v in df["requirement3_final_integrated_decision"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nFinal top 2 guides per gene after integrated Requirement 3 evidence:\n")
        for _, r in top.iterrows():
            f.write(
                f"- {r['target_gene']} rank {int(r['requirement3_final_integrated_rank_within_gene'])}: "
                f"{r['guide_id']} | "
                f"hits={int(r.get('total_hits_req3', 0))} | "
                f"Req2_burden={r.get('db_refined_total_burden', '')} | "
                f"Req3_mobile_burden={r.get('requirement3_mobile_context_burden', '')} | "
                f"direct_plasmid={r.get('direct_plasmid_evidence_strength', '')} | "
                f"direct_SCCmec={r.get('direct_sccmec_evidence_strength', '')} | "
                f"direct_mobile={r.get('direct_mobile_element_evidence_strength', '')} | "
                f"integrated_score={int(r['requirement3_integrated_evidence_score'])} | "
                f"final={r['requirement3_final_integrated_decision']}\n"
            )

        f.write("\nRequirement 3 status:\n")
        f.write(
            "Requirement 3 is complete as a direct-evidence-enhanced computational "
            "mapping layer. The analysis now includes plasmid/chromosome contig "
            "classification, SCCmec-like marker evidence for mecA neighborhoods, "
            "integron/insertion-sequence/transposon marker evidence, resistance-island-like "
            "context scoring, operon-proxy scoring, all-20 guide restoration, and zero-hit "
            "guide retention. Dedicated external typing tools such as PlasmidFinder, "
            "SCCmecFinder, IntegronFinder, ISfinder/ISEScan, or MobileElementFinder can "
            "still be described as optional future validation rather than required core "
            "pipeline steps.\n"
        )

    print(f"Wrote: {OUT_ALL20.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_TOP.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_DEPRIORITIZED.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nFinal top 2 guides per gene:")
    show_cols = [
        "target_gene",
        "requirement3_final_integrated_rank_within_gene",
        "guide_id",
        "total_hits_req3",
        "db_refined_total_burden",
        "requirement3_mobile_context_burden",
        "direct_plasmid_evidence_strength",
        "direct_sccmec_evidence_strength",
        "direct_mobile_element_evidence_strength",
        "requirement3_integrated_evidence_score",
        "requirement3_final_integrated_decision",
    ]
    show_cols = [c for c in show_cols if c in top.columns]
    print(top[show_cols].to_string(index=False))

    print("\nRequirement 3 final integrated confidence table completed.")


if __name__ == "__main__":
    main()
