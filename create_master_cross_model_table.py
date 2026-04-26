import pandas as pd

input_file = "results/cross_model_reference_ranked_guides.csv"
output_file = "results/master_cross_model_table.csv"

df = pd.read_csv(input_file)

master_df = pd.DataFrame({
    "gene": df["gene"],
    "position": df["position"],
    "guide": df["spacer"],
    "RS3_score": df["on_target_score"],
    "CFD_score": df["specificity_score"],
    "conservation": df["conservation_score"],
    "final_score": df["final_score_new"],
    "rank_RS3": df["rank_on_target"],
    "rank_CFD": df["rank_off_target"],
    "rank_final": df["rank_final"]
})

master_df = master_df.sort_values("rank_final")

master_df.to_csv(output_file, index=False)

print("✅ Clean master table created:", output_file)
print("Rows:", len(master_df))
print(master_df.head(10))
