import pandas as pd

df = pd.read_csv("results/final_comparison_table.csv")

print("Available columns:")
print(df.columns.tolist())

def find_col(possible_names):
    for name in possible_names:
        if name in df.columns:
            return name
    return None

guide_col = find_col(["guide_sequence", "spacer", "Target sequence", "target_sequence", "sequence"])
final_col = find_col(["final_score", "composite_score", "ensemble_score"])
mit_col = find_col(["mit_score", "MIT_score", "mit_guide_score", "MIT Guide Score", "mit_offtarget_score"])
crisot_col = find_col(["crisot_score", "CRISOT-Score", "CRISOT_score"])
cfd_col = find_col(["cfd_score", "CFD_score", "cfd_specificity_score"])
cons_col = find_col(["conservation_score", "panstrain_score", "pan_strain_score", "strain_conservation_score"])

needed = {
    "guide_col": guide_col,
    "final_col": final_col,
    "mit_col": mit_col,
    "crisot_col": crisot_col,
}

missing = [k for k, v in needed.items() if v is None]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

top_composite = df.sort_values(final_col, ascending=False).head(10)
top_mit = df.sort_values(mit_col, ascending=False).head(10)
top_crisot = df.sort_values(crisot_col, ascending=False).head(10)

comp_set = set(top_composite[guide_col])
mit_set = set(top_mit[guide_col])
crisot_set = set(top_crisot[guide_col])

print("\nTop 10 Overlap:")
print("Composite vs MIT:", len(comp_set & mit_set))
print("Composite vs CRISOT:", len(comp_set & crisot_set))

unique = comp_set - (mit_set | crisot_set)

print("\nGuides ONLY selected by your composite:")
for g in unique:
    print(g)

show_cols = [guide_col, final_col, mit_col, crisot_col]
for optional_col in [cfd_col, cons_col]:
    if optional_col is not None:
        show_cols.append(optional_col)

print("\nDetails of unique guides:")
print(df[df[guide_col].isin(unique)][show_cols])
