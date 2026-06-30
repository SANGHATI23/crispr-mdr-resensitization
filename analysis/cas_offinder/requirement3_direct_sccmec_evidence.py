#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HIT_CONTEXT = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/requirement3_hit_level_mobile_context_mapping.csv"

ALL20_DIRECT_PLASMID = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/final_requirement3_direct_evidence/requirement3_all20_with_direct_plasmid_evidence.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_sccmec_mapping"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_HIT = OUT_DIR / "requirement3_direct_sccmec_hit_level_mapping.csv"
OUT_GUIDE = OUT_DIR / "requirement3_direct_sccmec_guide_level_summary.csv"
OUT_ALL20 = OUT_DIR / "requirement3_all20_with_direct_plasmid_and_sccmec_evidence.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_direct_sccmec_mapping_summary.txt"


SCCMEC_MARKERS = [
    "sccmec",
    "staphylococcal cassette chromosome",
    "cassette chromosome",
    "chromosomal cassette",
    "cassette chromosome recombinase",
    "meca",
    "mecr1",
    "meci",
    "ccra",
    "ccrb",
    "ccrc",
    "orfx",
]


def contains_marker(text):
    text = str(text).lower()
    return any(marker in text for marker in SCCMEC_MARKERS)


def classify_strength(row):
    marker_hits = int(row.get("direct_sccmec_marker_hits", 0))
    proxy_hits = int(row.get("direct_sccmec_proxy_context_hits", 0))
    relevant_hits = int(row.get("direct_sccmec_relevant_hits", 0))

    if marker_hits >= 10 or proxy_hits >= 10 or relevant_hits >= 10:
        return "Strong_direct_SCCmec_like_evidence"
    if marker_hits >= 3 or proxy_hits >= 3 or relevant_hits >= 3:
        return "Moderate_direct_SCCmec_like_evidence"
    if marker_hits > 0 or proxy_hits > 0 or relevant_hits > 0:
        return "Low_direct_SCCmec_like_evidence"
    return "No_direct_SCCmec_like_evidence"


def main():
    print("Requirement 3 direct SCCmec evidence mapping started")

    if not HIT_CONTEXT.exists():
        raise FileNotFoundError(f"Missing hit context table: {HIT_CONTEXT}")

    if not ALL20_DIRECT_PLASMID.exists():
        raise FileNotFoundError(f"Missing all-20 direct plasmid table: {ALL20_DIRECT_PLASMID}")

    hits = pd.read_csv(HIT_CONTEXT)
    all20 = pd.read_csv(ALL20_DIRECT_PLASMID)

    print(f"Loaded hit context table: {hits.shape}")
    print(f"Loaded all-20 direct plasmid table: {all20.shape}")

    required = ["target_gene", "guide_id", "nearby_gene_products"]
    for col in required:
        if col not in hits.columns:
            raise ValueError(f"Missing required column in hit table: {col}")

    if "sccmec_context_flag" not in hits.columns:
        hits["sccmec_context_flag"] = 0

    if "nearby_sccmec_feature_count" not in hits.columns:
        hits["nearby_sccmec_feature_count"] = 0

    hits["direct_sccmec_marker_flag"] = hits["nearby_gene_products"].apply(
        lambda x: int(contains_marker(x))
    )

    hits["is_meca_target"] = hits["target_gene"].astype(str).str.lower().eq("meca").astype(int)

    hits["direct_sccmec_relevant_hit"] = (
        (hits["is_meca_target"] == 1)
        & (
            (hits["direct_sccmec_marker_flag"] > 0)
            | (hits["sccmec_context_flag"] > 0)
            | (hits["nearby_sccmec_feature_count"] > 0)
        )
    ).astype(int)

    hits.to_csv(OUT_HIT, index=False)

    guide = (
        hits.groupby(["target_gene", "guide_id"], dropna=False)
        .agg(
            total_hit_rows_sccmec=("guide_id", "count"),
            direct_sccmec_marker_hits=("direct_sccmec_marker_flag", "sum"),
            direct_sccmec_proxy_context_hits=("sccmec_context_flag", "sum"),
            direct_nearby_sccmec_feature_hits=("nearby_sccmec_feature_count", lambda x: int((x > 0).sum())),
            direct_sccmec_relevant_hits=("direct_sccmec_relevant_hit", "sum"),
        )
        .reset_index()
    )

    guide["direct_sccmec_evidence_strength"] = guide.apply(classify_strength, axis=1)

    guide.to_csv(OUT_GUIDE, index=False)

    keep_cols = [
        "target_gene",
        "guide_id",
        "direct_sccmec_marker_hits",
        "direct_sccmec_proxy_context_hits",
        "direct_nearby_sccmec_feature_hits",
        "direct_sccmec_relevant_hits",
        "direct_sccmec_evidence_strength",
    ]

    merged = all20.merge(
        guide[keep_cols],
        on=["target_gene", "guide_id"],
        how="left",
    )

    zero_cols = [
        "direct_sccmec_marker_hits",
        "direct_sccmec_proxy_context_hits",
        "direct_nearby_sccmec_feature_hits",
        "direct_sccmec_relevant_hits",
    ]

    for col in zero_cols:
        if col not in merged.columns:
            merged[col] = 0
        merged[col] = merged[col].fillna(0).astype(int)

    merged["direct_sccmec_evidence_strength"] = merged[
        "direct_sccmec_evidence_strength"
    ].fillna("Zero_hit_no_direct_SCCmec_like_evidence")

    merged.to_csv(OUT_ALL20, index=False)

    mecA_guide = guide[guide["target_gene"].astype(str).str.lower() == "meca"].copy()

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Direct SCCmec Evidence Summary\n")
        f.write("===========================================\n\n")

        f.write(f"Input hit context table: {HIT_CONTEXT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Input all-20 direct plasmid table: {ALL20_DIRECT_PLASMID.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Hit rows evaluated: {len(hits)}\n")
        f.write(f"Guide rows summarized before all-20 merge: {len(guide)}\n")
        f.write(f"Final all-20 rows after merge: {len(merged)}\n")
        f.write(f"Unique guides after merge: {merged['guide_id'].nunique()}\n\n")

        f.write("Overall SCCmec evidence strength counts:\n")
        for k, v in merged["direct_sccmec_evidence_strength"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nmecA guide-level SCCmec evidence:\n")
        if mecA_guide.empty:
            f.write("- No mecA guide rows found.\n")
        else:
            for _, r in mecA_guide.iterrows():
                f.write(
                    f"- {r['guide_id']} | "
                    f"total_hits={int(r['total_hit_rows_sccmec'])} | "
                    f"sccmec_marker_hits={int(r['direct_sccmec_marker_hits'])} | "
                    f"sccmec_proxy_context_hits={int(r['direct_sccmec_proxy_context_hits'])} | "
                    f"nearby_sccmec_feature_hits={int(r['direct_nearby_sccmec_feature_hits'])} | "
                    f"relevant_mecA_sccmec_hits={int(r['direct_sccmec_relevant_hits'])} | "
                    f"strength={r['direct_sccmec_evidence_strength']}\n"
                )

        f.write("\nOutput files:\n")
        f.write(f"- {OUT_HIT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_GUIDE.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_ALL20.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_SUMMARY.relative_to(PROJECT_ROOT)}\n")

    print(f"Wrote: {OUT_HIT.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_GUIDE.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_ALL20.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nmecA SCCmec evidence:")
    if mecA_guide.empty:
        print("No mecA guide rows found.")
    else:
        print(mecA_guide.to_string(index=False))

    print("\nRequirement 3 direct SCCmec evidence mapping completed.")


if __name__ == "__main__":
    main()
