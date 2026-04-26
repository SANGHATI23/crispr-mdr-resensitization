import pandas as pd

master_file = "results/master_cross_model_table.csv"
mit_file = "results_external/mit/mit_guide_scores.csv"

output_file = "results/master_cross_model_with_mit.csv"

master = pd.read_csv(master_file)
mit = pd.read_csv(mit_file)

merged = master.merge(mit, on="guide", how="left")

merged.to_csv(output_file, index=False)

print("✅ Final table with MIT:", output_file)
print(merged.head())
