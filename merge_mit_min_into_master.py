import pandas as pd

master_file = "results/master_cross_model_with_mit.csv"
mit_min_file = "results_external/mit/mit_guide_scores_min.csv"
output_file = "results/master_cross_model_with_mit_min.csv"

master = pd.read_csv(master_file)
mit_min = pd.read_csv(mit_min_file)

merged = master.merge(mit_min, on="guide", how="left")

merged.to_csv(output_file, index=False)

print("✅ Final table with mean + worst-case MIT:", output_file)
print("Rows:", len(merged))
print(merged.head())
