import pandas as pd

comp = pd.read_csv("results/master_cross_model_with_mit.csv")
crisot = pd.read_csv("results_external/crisot/crisot_rescore_output.csv")

print("MASTER COLUMNS:", comp.columns.tolist())
print("CRISOT COLUMNS:", crisot.columns.tolist())

# Find guide/spacer column in master file
possible_guide_cols = [
    "spacer",
    "guide",
    "guide_sequence",
    "spacer_sequence",
    "sequence",
    "Target sequence",
    "target_sequence",
]

master_guide_col = None
for c in possible_guide_cols:
    if c in comp.columns:
        master_guide_col = c
        break

if master_guide_col is None:
    raise ValueError("Could not find guide/spacer column in master file.")

# CRISOT guide column
crisot_guide_col = "Target sequence"

if crisot_guide_col not in crisot.columns:
    raise ValueError("Could not find 'Target sequence' column in CRISOT file.")

# Rename for merge
comp = comp.rename(columns={master_guide_col: "guide_sequence"})
crisot = crisot.rename(columns={
    "Target sequence": "guide_sequence",
    "CRISOT-Score": "crisot_score",
    "CRISOT-Spec": "crisot_spec",
})

# Clean sequence text
comp["guide_sequence"] = comp["guide_sequence"].astype(str).str.upper().str.strip()
crisot["guide_sequence"] = crisot["guide_sequence"].astype(str).str.upper().str.strip()

# Merge
df = comp.merge(
    crisot[["guide_sequence", "crisot_score", "crisot_spec"]],
    on="guide_sequence",
    how="left"
)

# Save full merged version first
df.to_csv("results/final_comparison_table.csv", index=False)

print("Saved: results/final_comparison_table.csv")
print("Rows:", len(df))
print("CRISOT matched rows:", df["crisot_score"].notna().sum())
print(df.head())
