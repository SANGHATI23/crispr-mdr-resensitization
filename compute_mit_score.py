import pandas as pd
import numpy as np

# -----------------------------
# Load pair-level CFD file
# -----------------------------
input_file = "results_external/cfd/cfd_scored.csv"
output_file = "results_external/mit/mit_scored_pairs.csv"

df = pd.read_csv(input_file)

# -----------------------------
# MIT mismatch position weights
# (standard approximate weights)
# -----------------------------
position_weights = {
    1: 0.0, 2: 0.0, 3: 0.014, 4: 0.0, 5: 0.0,
    6: 0.395, 7: 0.317, 8: 0.0, 9: 0.389, 10: 0.079,
    11: 0.445, 12: 0.508, 13: 0.613, 14: 0.851,
    15: 0.732, 16: 0.828, 17: 0.615, 18: 0.804,
    19: 0.685, 20: 0.583
}

def compute_mit(spacer, target):
    spacer = spacer[:20]
    target = target[:20]

    mismatches = []
    
    for i in range(20):
        if spacer[i] != target[i]:
            mismatches.append(i+1)

    if len(mismatches) == 0:
        return 1.0

    score = 1.0

    for pos in mismatches:
        w = position_weights.get(pos, 0.5)
        score *= (1 - w)

    # mismatch count penalty
    score *= 1 / (len(mismatches) ** 2)

    return score

# -----------------------------
# Compute MIT score
# -----------------------------
df["MIT_pair_score"] = df.apply(
    lambda row: compute_mit(row["spacer"], row["off_target"]),
    axis=1
)

# -----------------------------
# Save pair-level output
# -----------------------------
import os
os.makedirs("results_external/mit", exist_ok=True)

df.to_csv(output_file, index=False)

print("✅ MIT pair scores created:", output_file)
print("Rows:", len(df))
