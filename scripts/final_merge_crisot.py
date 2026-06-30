import pandas as pd

comp = pd.read_csv("results/master_cross_model_with_mit.csv")
crisot = pd.read_csv("results_external/crisot/crisot_master_rescore_output.csv")

# Detect guide column
possible_cols = ["guide", "guide_sequence", "spacer", "sequence"]
guide_col = next((c for c in possible_cols if c in comp.columns), None)

comp = comp.rename(columns={guide_col: "guide_sequence"})
comp["guide_sequence"] = comp["guide_sequence"].str.upper().str.strip()

# CRISOT → remove PAM
crisot["guide_sequence"] = crisot["Target sequence"].str[:20].str.upper().str.strip()
crisot = crisot.rename(columns={"CRISOT-Score": "crisot_score"})

df = comp.merge(
    crisot[["guide_sequence", "crisot_score"]],
    on="guide_sequence",
    how="left"
)

df.to_csv("results/final_comparison_table.csv", index=False)

print("Merged successfully")
print("Matched CRISOT rows:", df["crisot_score"].notna().sum())
print(df.head())
