import pandas as pd

comp = pd.read_csv("results/master_cross_model_with_mit.csv")
crisot = pd.read_csv("results_external/crisot/crisot_rescore_output.csv")

print("MASTER COLUMNS:")
print(comp.columns.tolist())

possible_cols = [
    "guide_sequence",
    "spacer",
    "guide",
    "sequence",
    "target",
    "Target sequence",
    "target_sequence",
    "sgRNA",
    "sgRNA_sequence",
    "gRNA",
    "gRNA_sequence"
]

guide_col = None
for c in possible_cols:
    if c in comp.columns:
        guide_col = c
        break

if guide_col is None:
    raise ValueError("Could not find guide column. Look at MASTER COLUMNS above and tell me which column contains the 20nt guide.")

print("Using guide column:", guide_col)

comp = comp.rename(columns={guide_col: "guide_sequence"})

comp["guide_sequence"] = (
    comp["guide_sequence"]
    .astype(str)
    .str.upper()
    .str.strip()
)

crisot["guide_sequence"] = (
    crisot["Target sequence"]
    .astype(str)
    .str.upper()
    .str.strip()
    .str[:20]
)

crisot = crisot.rename(columns={"CRISOT-Score": "crisot_score"})

df = comp.merge(
    crisot[["guide_sequence", "crisot_score"]],
    on="guide_sequence",
    how="left"
)

df.to_csv("results/final_comparison_table.csv", index=False)

print("Fixed CRISOT merge done")
print("Rows:", len(df))
print("Matched CRISOT rows:", df["crisot_score"].notna().sum())
print(df[["guide_sequence", "crisot_score"]].head(10))
