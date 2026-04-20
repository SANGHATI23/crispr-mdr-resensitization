import pandas as pd

# Read the three main Phase 2 files
w2 = pd.read_csv("weight_sensitivity_results/phase2_weight_summary_all_variations.csv")
g2 = pd.read_csv("weight_sensitivity_results/phase2_per_gene_summary_all_variations.csv")
a2 = pd.read_csv("weight_sensitivity_results/phase2_all_weighted_guides_all_variations.csv")

print("\n=== PHASE 2 WEIGHT SUMMARY ===")
print(w2.head())
print("\nColumns:")
print(w2.columns.tolist())

print("\n=== PHASE 2 PER GENE SUMMARY ===")
print(g2.head())
print("\nColumns:")
print(g2.columns.tolist())

print("\n=== PHASE 2 ALL GUIDES ===")
print(a2.head())
print("\nColumns:")
print(a2.columns.tolist())