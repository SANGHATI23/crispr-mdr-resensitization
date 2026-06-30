import pandas as pd
from pathlib import Path

ROOT = Path(".")
OUT = ROOT / "results_bibm_strengthening" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

INPUT = ROOT / "results_mlcb" / "tables" / "functional_offtarget_risk_ranked_guides.csv"

if not INPUT.exists():
    raise FileNotFoundError(f"Missing input file: {INPUT}")

df = pd.read_csv(INPUT)

print("Loaded:", INPUT)
print("Shape:", df.shape)
print("Columns:", list(df.columns))

required = [
    "mlcb_gene",
    "mlcb_spacer",
    "mlcb_base_score_scaled",
    "crisot_wt_specificity",
    "annotation_derived_functional_risk_score",
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Standardized columns
df["gene"] = df["mlcb_gene"].astype(str)
df["guide_id"] = (
    df["mlcb_gene"].astype(str)
    + "_"
    + df["mlcb_spacer"].astype(str).str[:10]
    + "_row"
    + df.index.astype(str)
)

# Component scores
# activity_score: already scaled 0-1
# specificity_score: CRISOT WT specificity 0-1
# functional_safety_score: inverted functional risk, so higher = safer
df["activity_score"] = pd.to_numeric(df["mlcb_base_score_scaled"], errors="coerce")
df["specificity_score"] = pd.to_numeric(df["crisot_wt_specificity"], errors="coerce")
df["functional_risk_score"] = pd.to_numeric(
    df["annotation_derived_functional_risk_score"],
    errors="coerce"
)
df["functional_safety_score"] = 1 - df["functional_risk_score"]

df = df.dropna(
    subset=[
        "gene",
        "guide_id",
        "activity_score",
        "specificity_score",
        "functional_safety_score",
    ]
).copy()

if df.empty:
    raise ValueError("No usable rows remained after numeric conversion.")

# Clip only for safety
for c in ["activity_score", "specificity_score", "functional_safety_score"]:
    df[c] = df[c].clip(0, 1)

extreme_weights = [
    (80, 10, 10),
    (10, 80, 10),
    (10, 10, 80),
    (70, 20, 10),
    (70, 10, 20),
    (20, 70, 10),
    (10, 70, 20),
    (20, 10, 70),
    (10, 20, 70),
    (90, 5, 5),
    (5, 90, 5),
    (5, 5, 90),
]

# Reference FOTR-v2 balanced score
df["score_40_30_30"] = (
    0.40 * df["activity_score"]
    + 0.30 * df["specificity_score"]
    + 0.30 * df["functional_safety_score"]
)

ref_top = (
    df.sort_values(["gene", "score_40_30_30"], ascending=[True, False])
      .groupby("gene")
      .head(1)[["gene", "guide_id", "mlcb_spacer"]]
      .rename(columns={
          "guide_id": "reference_top1",
          "mlcb_spacer": "reference_top1_spacer",
      })
)

ref_top20 = set(
    df.sort_values("score_40_30_30", ascending=False)
      .head(20)["guide_id"]
      .astype(str)
)

rows = []
detail_rows = []

for aw, sw, fw in extreme_weights:
    label = f"{aw}_{sw}_{fw}"
    score_col = f"score_{label}"

    df[score_col] = (
        (aw / 100) * df["activity_score"]
        + (sw / 100) * df["specificity_score"]
        + (fw / 100) * df["functional_safety_score"]
    )

    top = (
        df.sort_values(["gene", score_col], ascending=[True, False])
          .groupby("gene")
          .head(1)[["gene", "guide_id", "mlcb_spacer"]]
          .rename(columns={
              "guide_id": "extreme_top1",
              "mlcb_spacer": "extreme_top1_spacer",
          })
    )

    merged = ref_top.merge(top, on="gene", how="left")

    merged["weight_configuration"] = label
    merged["top1_stable"] = (
        merged["reference_top1"].astype(str)
        == merged["extreme_top1"].astype(str)
    )

    stable_genes = merged["top1_stable"].sum()
    total_genes = merged.shape[0]

    top20 = set(
        df.sort_values(score_col, ascending=False)
          .head(20)["guide_id"]
          .astype(str)
    )

    union = ref_top20 | top20
    jaccard = len(ref_top20 & top20) / len(union) if union else 0

    rows.append({
        "source_file": str(INPUT),
        "weight_configuration": label,
        "activity_weight": aw,
        "specificity_weight": sw,
        "functional_safety_weight": fw,
        "stable_gene_top1_count": int(stable_genes),
        "total_genes": int(total_genes),
        "all_gene_top1_stable": bool(stable_genes == total_genes),
        "top20_jaccard_vs_40_30_30": round(jaccard, 4),
    })

    detail_rows.append(merged)

summary = pd.DataFrame(rows)
details = pd.concat(detail_rows, ignore_index=True)

summary_path = OUT / "extreme_weight_robustness.csv"
details_path = OUT / "extreme_weight_top1_details.csv"
scored_path = OUT / "extreme_weight_scored_guides.csv"

summary.to_csv(summary_path, index=False)
details.to_csv(details_path, index=False)
df.to_csv(scored_path, index=False)

print("\nExtreme robustness summary:")
print(summary.to_string(index=False))

print("\nTop-1 details:")
print(details.to_string(index=False))

print(f"\nWrote: {summary_path}")
print(f"Wrote: {details_path}")
print(f"Wrote: {scored_path}")
