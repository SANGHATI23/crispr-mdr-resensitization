import pandas as pd
from pathlib import Path

df = pd.read_csv("results/master_cross_model_with_mit.csv")

print("Available columns:")
print(df.columns.tolist())

possible_cols = [
    "guide_sequence",
    "spacer",
    "guide",
    "sequence",
    "target",
    "Target sequence",
    "target_sequence",
    "sgRNA",
    "sgRNA_sequence",
    "gRNA",
    "gRNA_sequence",
]

guide_col = None
for c in possible_cols:
    if c in df.columns:
        guide_col = c
        break

if guide_col is None:
    raise ValueError("Could not auto-detect guide column. Check the printed columns above.")

print("Using guide column:", guide_col)

outdir = Path("results_external/crisot")
outdir.mkdir(parents=True, exist_ok=True)

guides = (
    df[guide_col]
    .astype(str)
    .str.upper()
    .str.strip()
    .dropna()
    .drop_duplicates()
)

targets_23nt = guides + "TGG"

guides.to_csv(outdir / "crisot_master_guides_20nt.txt", index=False, header=False)
targets_23nt.to_csv(outdir / "crisot_master_targets_23nt.txt", index=False, header=False)

print("Saved:")
print(outdir / "crisot_master_guides_20nt.txt")
print(outdir / "crisot_master_targets_23nt.txt")
print("Unique guides:", len(guides))
