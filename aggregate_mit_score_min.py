import pandas as pd

input_file = "results_external/mit/mit_scored_pairs.csv"
output_file = "results_external/mit/mit_guide_scores_min.csv"

df = pd.read_csv(input_file)

guide_df = df.groupby("spacer")["MIT_pair_score"].min().reset_index()

guide_df = guide_df.rename(columns={
    "spacer": "guide",
    "MIT_pair_score": "MIT_score_min"
})

guide_df.to_csv(output_file, index=False)

print("✅ Created:", output_file)
print("Rows:", len(guide_df))
print(guide_df.head())
