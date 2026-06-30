from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[3]

INPUT_GUIDES = ROOT / "results_journal_extension" / "unseen_genes" / "unseen_gene_candidate_guides_all.csv"

OUT_DIR = ROOT / "results_journal_extension" / "unseen_genes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "unseen_gene_candidate_guides_ranked_v1.csv"
OUT_SUMMARY = OUT_DIR / "unseen_gene_candidate_guides_ranked_v1_summary.txt"


def region_weight(region):
    if region == "early_coding_region":
        return 1.00
    if region == "middle_coding_region":
        return 0.85
    if region == "late_coding_region":
        return 0.70
    return 0.75


def gc_weight(gc_class):
    if gc_class == "Acceptable_GC":
        return 1.00
    if gc_class in ["Low_GC", "High_GC"]:
        return 0.70
    return 0.80


def classify_priority(score, gc_class):
    if gc_class != "Acceptable_GC":
        return "Lower_priority_gc_extreme"
    if score >= 85:
        return "High_priority_unseen_candidate"
    if score >= 70:
        return "Moderate_priority_unseen_candidate"
    return "Lower_priority_unseen_candidate"


def main():
    df = pd.read_csv(INPUT_GUIDES)

    if df.empty:
        raise ValueError("No unseen candidate guides found. Run enumerate_unseen_gene_guides.py first.")

    df["region_weight_v1"] = df["target_region_bin"].apply(region_weight)
    df["gc_weight_v1"] = df["gc_class"].apply(gc_weight)

    # First-pass prospective unseen-gene score.
    # This is intentionally conservative and transparent.
    df["unseen_gene_priority_score_v1"] = (
        0.70 * df["simple_on_target_proxy"] +
        20.0 * df["region_weight_v1"] +
        10.0 * df["gc_weight_v1"]
    )

    # Normalize to 0-100-ish scale
    max_possible = 0.70 * 100 + 20.0 * 1.0 + 10.0 * 1.0
    df["unseen_gene_priority_score_v1"] = (df["unseen_gene_priority_score_v1"] / max_possible) * 100
    df["unseen_gene_priority_score_v1"] = df["unseen_gene_priority_score_v1"].round(2)

    df["unseen_gene_rank_v1"] = df.groupby("gene")["unseen_gene_priority_score_v1"].rank(
        method="min", ascending=False
    ).astype(int)

    df["unseen_gene_recommendation_v1"] = df.apply(
        lambda r: classify_priority(r["unseen_gene_priority_score_v1"], r["gc_class"]),
        axis=1
    )

    df = df.sort_values(["gene", "unseen_gene_rank_v1", "position_1based"]).reset_index(drop=True)

    df.to_csv(OUT_TABLE, index=False)

    lines = []
    lines.append("Unseen Gene Candidate Guide Ranking v1 Summary")
    lines.append("=============================================")
    lines.append("")
    lines.append(f"Input guides: {INPUT_GUIDES}")
    lines.append(f"Output ranked table: {OUT_TABLE}")
    lines.append(f"Total ranked guides: {len(df)}")
    lines.append("")

    lines.append("Recommendation counts:")
    for k, v in df["unseen_gene_recommendation_v1"].value_counts().items():
        lines.append(f"- {k}: {v}")

    lines.append("")
    lines.append("Top 10 guides by gene:")
    for gene, gdf in df.groupby("gene"):
        lines.append(f"{gene}:")
        top = gdf.sort_values("unseen_gene_rank_v1").head(10)
        for _, r in top.iterrows():
            lines.append(
                f"  - rank={r['unseen_gene_rank_v1']} | spacer={r['spacer']} | PAM={r['pam']} | "
                f"strand={r['strand']} | pos={r['position_1based']} | "
                f"region={r['target_region_bin']} | GC={r['gc_content']} | "
                f"score={r['unseen_gene_priority_score_v1']} | "
                f"label={r['unseen_gene_recommendation_v1']}"
            )

    lines.append("")
    lines.append("Interpretation:")
    lines.append(
        "This v1 ranking provides a transparent first-pass prioritization of unseen AMR-gene guides. "
        "It uses sequence suitability, GC status, and target-region preference only. It should be treated as "
        "a prospective-screening layer, not as final FOTR v2-equivalent scoring until RNA accessibility, external "
        "off-target screening, conservation, and functional annotation are added."
    )

    OUT_SUMMARY.write_text("\n".join(lines))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
