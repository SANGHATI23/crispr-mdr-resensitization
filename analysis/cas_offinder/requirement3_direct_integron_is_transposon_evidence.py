#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HIT_CONTEXT = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/requirement3_hit_level_mobile_context_mapping.csv"

ALL20_WITH_SCCMEC = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_sccmec_mapping/requirement3_all20_with_direct_plasmid_and_sccmec_evidence.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_integron_is_transposon_mapping"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_HIT = OUT_DIR / "requirement3_direct_integron_is_transposon_hit_level_mapping.csv"
OUT_GUIDE = OUT_DIR / "requirement3_direct_integron_is_transposon_guide_level_summary.csv"
OUT_ALL20 = OUT_DIR / "requirement3_all20_with_direct_plasmid_sccmec_integron_is_transposon_evidence.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_direct_integron_is_transposon_mapping_summary.txt"


INTEGRON_MARKERS = [
    "integron",
    "inti1",
    "inti2",
    "inti3",
    "integrase inti",
    "gene cassette",
    "attc",
    "qac",
    "sul1",
]

INSERTION_SEQUENCE_MARKERS = [
    "insertion sequence",
    "insertion element",
    "is family",
    "is3",
    "is5",
    "is6",
    "is30",
    "is110",
    "is256",
    "is630",
    "is982",
    "is1182",
    "is1380",
    "is481",
    "is4 family",
    "is5 family",
    "is6 family",
    "transposase",
]

TRANSPOSON_MARKERS = [
    "transposon",
    "tn3",
    "tn7",
    "tn4401",
    "tn1546",
    "tn916",
    "tn21",
    "tn10",
    "conjugative transposon",
    "resolvase",
    "recombinase",
    "invertase",
    "transposition",
]

CONJUGATION_MARKERS = [
    "conjugative",
    "conjugal",
    "tra protein",
    "transfer protein",
    "type iv secretion",
    "relaxase",
    "mobilization",
    "mobilisation",
    "mob",
]


def contains_any(text, markers):
    text = str(text).lower()
    return any(marker in text for marker in markers)


def classify_strength(row):
    integron_hits = int(row.get("direct_integron_marker_hits", 0))
    is_hits = int(row.get("direct_insertion_sequence_marker_hits", 0))
    transposon_hits = int(row.get("direct_transposon_marker_hits", 0))
    conjugation_hits = int(row.get("direct_conjugation_marker_hits", 0))
    total = integron_hits + is_hits + transposon_hits + conjugation_hits

    if total >= 50:
        return "Strong_direct_mobile_element_evidence"
    if total >= 10:
        return "Moderate_direct_mobile_element_evidence"
    if total > 0:
        return "Low_direct_mobile_element_evidence"
    return "No_direct_mobile_element_evidence"


def main():
    print("Requirement 3 direct integron/IS/transposon evidence mapping started")

    if not HIT_CONTEXT.exists():
        raise FileNotFoundError(f"Missing hit context table: {HIT_CONTEXT}")

    if not ALL20_WITH_SCCMEC.exists():
        raise FileNotFoundError(f"Missing all-20 plasmid+SCCmec table: {ALL20_WITH_SCCMEC}")

    hits = pd.read_csv(HIT_CONTEXT)
    all20 = pd.read_csv(ALL20_WITH_SCCMEC)

    print(f"Loaded hit context table: {hits.shape}")
    print(f"Loaded all-20 plasmid+SCCmec table: {all20.shape}")

    required = ["target_gene", "guide_id", "nearby_gene_products"]
    for col in required:
        if col not in hits.columns:
            raise ValueError(f"Missing required column in hit table: {col}")

    feature_text = hits["nearby_gene_products"].fillna("").astype(str)

    hits["direct_integron_marker_flag"] = feature_text.apply(
        lambda x: int(contains_any(x, INTEGRON_MARKERS))
    )
    hits["direct_insertion_sequence_marker_flag"] = feature_text.apply(
        lambda x: int(contains_any(x, INSERTION_SEQUENCE_MARKERS))
    )
    hits["direct_transposon_marker_flag"] = feature_text.apply(
        lambda x: int(contains_any(x, TRANSPOSON_MARKERS))
    )
    hits["direct_conjugation_marker_flag"] = feature_text.apply(
        lambda x: int(contains_any(x, CONJUGATION_MARKERS))
    )

    hits["direct_mobile_element_any_flag"] = (
        (
            hits["direct_integron_marker_flag"]
            + hits["direct_insertion_sequence_marker_flag"]
            + hits["direct_transposon_marker_flag"]
            + hits["direct_conjugation_marker_flag"]
        ) > 0
    ).astype(int)

    hits.to_csv(OUT_HIT, index=False)

    guide = (
        hits.groupby(["target_gene", "guide_id"], dropna=False)
        .agg(
            total_hit_rows_mobile=("guide_id", "count"),
            direct_integron_marker_hits=("direct_integron_marker_flag", "sum"),
            direct_insertion_sequence_marker_hits=("direct_insertion_sequence_marker_flag", "sum"),
            direct_transposon_marker_hits=("direct_transposon_marker_flag", "sum"),
            direct_conjugation_marker_hits=("direct_conjugation_marker_flag", "sum"),
            direct_any_mobile_element_marker_hits=("direct_mobile_element_any_flag", "sum"),
        )
        .reset_index()
    )

    guide["direct_mobile_element_evidence_strength"] = guide.apply(classify_strength, axis=1)

    guide.to_csv(OUT_GUIDE, index=False)

    keep_cols = [
        "target_gene",
        "guide_id",
        "direct_integron_marker_hits",
        "direct_insertion_sequence_marker_hits",
        "direct_transposon_marker_hits",
        "direct_conjugation_marker_hits",
        "direct_any_mobile_element_marker_hits",
        "direct_mobile_element_evidence_strength",
    ]

    merged = all20.merge(
        guide[keep_cols],
        on=["target_gene", "guide_id"],
        how="left",
    )

    zero_cols = [
        "direct_integron_marker_hits",
        "direct_insertion_sequence_marker_hits",
        "direct_transposon_marker_hits",
        "direct_conjugation_marker_hits",
        "direct_any_mobile_element_marker_hits",
    ]

    for col in zero_cols:
        if col not in merged.columns:
            merged[col] = 0
        merged[col] = merged[col].fillna(0).astype(int)

    merged["direct_mobile_element_evidence_strength"] = merged[
        "direct_mobile_element_evidence_strength"
    ].fillna("Zero_hit_no_direct_mobile_element_evidence")

    merged.to_csv(OUT_ALL20, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Direct Integron / IS / Transposon Evidence Summary\n")
        f.write("================================================================\n\n")

        f.write(f"Input hit context table: {HIT_CONTEXT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Input all-20 plasmid+SCCmec table: {ALL20_WITH_SCCMEC.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Hit rows evaluated: {len(hits)}\n")
        f.write(f"Guide rows summarized before all-20 merge: {len(guide)}\n")
        f.write(f"Final all-20 rows after merge: {len(merged)}\n")
        f.write(f"Unique guides after merge: {merged['guide_id'].nunique()}\n\n")

        f.write("Hit-level direct mobile marker counts:\n")
        for col in [
            "direct_integron_marker_flag",
            "direct_insertion_sequence_marker_flag",
            "direct_transposon_marker_flag",
            "direct_conjugation_marker_flag",
            "direct_mobile_element_any_flag",
        ]:
            f.write(f"- {col}: {int(hits[col].sum())}\n")

        f.write("\nGuide-level direct mobile-element evidence strength counts:\n")
        for k, v in merged["direct_mobile_element_evidence_strength"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nGuides with any direct mobile-element marker evidence:\n")
        subset = merged[merged["direct_any_mobile_element_marker_hits"] > 0].copy()
        if subset.empty:
            f.write("- None\n")
        else:
            subset = subset.sort_values("direct_any_mobile_element_marker_hits", ascending=False)
            for _, r in subset.iterrows():
                f.write(
                    f"- {r['target_gene']} | {r['guide_id']} | "
                    f"any_mobile_hits={int(r['direct_any_mobile_element_marker_hits'])} | "
                    f"integron={int(r['direct_integron_marker_hits'])} | "
                    f"IS={int(r['direct_insertion_sequence_marker_hits'])} | "
                    f"transposon={int(r['direct_transposon_marker_hits'])} | "
                    f"conjugation={int(r['direct_conjugation_marker_hits'])} | "
                    f"strength={r['direct_mobile_element_evidence_strength']}\n"
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

    print("\nDirect mobile-element evidence strength counts:")
    print(merged["direct_mobile_element_evidence_strength"].value_counts().to_string())

    print("\nRequirement 3 direct integron/IS/transposon evidence mapping completed.")


if __name__ == "__main__":
    main()
