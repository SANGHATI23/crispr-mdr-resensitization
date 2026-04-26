import pandas as pd

input_file = "results_external/mit/mit_scored_pairs.csv"
output_file = "results_external/mit/mit_guide_scores.csv"

df = pd.read_csv(input_file)

# lower MIT = safer (less off-target cutting)
guide_df = df.groupby("spacer")["MIT_pair_score"].mean().reset_index()

guide_df = guide_df.rename(columns={
    "spacer": "guide",
    "MIT_pair_score": "MIT_score"
})

guide_df.to_csv(output_file, index=False)

print("✅ Guide-level MIT scores:", output_file)
print(guide_df.head())
