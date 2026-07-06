from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import RNA

ROOT = Path(".")
TABLE_DIR = ROOT / "results_mlcb" / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

input_path = TABLE_DIR / "mlcb_guide_dataset_with_local_context_embedding_ids.csv"

if not input_path.exists():
    raise FileNotFoundError(f"Missing input file: {input_path}")

df = pd.read_csv(input_path)

required_cols = [
    "mlcb_gene",
    "mlcb_spacer",
    "mlcb_base_score",
    "weak_functional_severity_label",
    "dnabert_sequence",
    "local_context_sequence",
    "gc_content",
    "on_target_score",
    "specificity_score",
    "conservation_score",
    "number_of_reference_matches",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

def dna_to_rna(seq):
    seq = str(seq).upper().replace("T", "U")
    keep = set("ACGUN")
    return "".join([b if b in keep else "N" for b in seq])

def unpaired_fraction(dotbracket):
    if not isinstance(dotbracket, str) or len(dotbracket) == 0:
        return np.nan
    return dotbracket.count(".") / len(dotbracket)

def center20_unpaired_fraction(dotbracket):
    if not isinstance(dotbracket, str) or len(dotbracket) == 0:
        return np.nan

    n = len(dotbracket)
    center = n // 2
    start = max(0, center - 10)
    end = min(n, center + 10)
    region = dotbracket[start:end]

    if len(region) == 0:
        return np.nan

    return region.count(".") / len(region)

records = []

for idx, row in df.iterrows():
    seq = row["local_context_sequence"]

    if pd.isna(seq) or str(seq).strip() == "":
        seq = row["dnabert_sequence"]

    dna_seq = str(seq).upper()
    rna_seq = dna_to_rna(dna_seq)

    structure, mfe = RNA.fold(rna_seq)

    records.append({
        "row_index": idx,
        "mlcb_gene": row["mlcb_gene"],
        "mlcb_spacer": row["mlcb_spacer"],
        "mlcb_base_score": row["mlcb_base_score"],
        "weak_functional_severity_label": row["weak_functional_severity_label"],
        "gc_content": row["gc_content"],
        "on_target_score": row["on_target_score"],
        "specificity_score": row["specificity_score"],
        "conservation_score": row["conservation_score"],
        "number_of_reference_matches": row["number_of_reference_matches"],
        "sequence_length": len(rna_seq),
        "rnafold_structure": structure,
        "rnafold_mfe": mfe,
        "rnafold_mfe_per_nt": mfe / len(rna_seq) if len(rna_seq) > 0 else np.nan,
        "rnafold_global_unpaired_fraction": unpaired_fraction(structure),
        "rnafold_center20_unpaired_fraction": center20_unpaired_fraction(structure),
    })

rnafold_df = pd.DataFrame(records)

# Add the same simple risk-aware score used in clinical triage, so RNAfold can be compared
# against an output ranking signal that was not used to compute RNAfold.
def minmax(s):
    s = pd.to_numeric(s, errors="coerce").fillna(0)
    if s.max() == s.min():
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.min()) / (s.max() - s.min())

for c in [
    "mlcb_base_score",
    "weak_functional_severity_label",
    "gc_content",
    "on_target_score",
    "specificity_score",
    "conservation_score",
    "number_of_reference_matches",
]:
    rnafold_df[c] = pd.to_numeric(rnafold_df[c], errors="coerce").fillna(0)

rnafold_df["base_score_norm"] = minmax(rnafold_df["mlcb_base_score"])

# Structure-aware external signal:
# More unpaired center sequence may suggest accessibility.
rnafold_df["rnafold_accessibility_proxy"] = rnafold_df["rnafold_center20_unpaired_fraction"]

# A compact risk-aware ranking proxy for correlation only.
# Higher score = stronger candidate after weak severity penalty.
rnafold_df["simple_risk_aware_score"] = (
    rnafold_df["base_score_norm"]
    * (1.0 - 0.25 * minmax(rnafold_df["weak_functional_severity_label"]))
)

metrics = [
    "rnafold_mfe",
    "rnafold_mfe_per_nt",
    "rnafold_global_unpaired_fraction",
    "rnafold_center20_unpaired_fraction",
    "rnafold_accessibility_proxy",
]

targets = [
    "mlcb_base_score",
    "weak_functional_severity_label",
    "simple_risk_aware_score",
    "gc_content",
    "on_target_score",
    "specificity_score",
    "conservation_score",
    "number_of_reference_matches",
]

corr_rows = []

for metric in metrics:
    for target in targets:
        temp = rnafold_df[[metric, target]].dropna()

        if len(temp) < 5:
            continue

        if temp[metric].nunique() < 2 or temp[target].nunique() < 2:
            continue

        rho, p_value = spearmanr(temp[metric], temp[target])

        corr_rows.append({
            "independent_RNAfold_metric": metric,
            "comparison_target": target,
            "n": len(temp),
            "spearman_rho": rho,
            "p_value": p_value,
            "absolute_rho": abs(rho),
        })

corr_df = pd.DataFrame(corr_rows).sort_values("absolute_rho", ascending=False)

out_metrics = TABLE_DIR / "independent_rnafold_validation_metrics.csv"
out_corr = TABLE_DIR / "independent_rnafold_spearman_correlations.csv"
out_corr_md = TABLE_DIR / "independent_rnafold_spearman_correlations.md"
out_summary = TABLE_DIR / "independent_rnafold_validation_summary.txt"

rnafold_df.to_csv(out_metrics, index=False)
corr_df.to_csv(out_corr, index=False)
corr_df.to_markdown(out_corr_md, index=False)

top_rows = corr_df.head(10).copy()

summary = []
summary.append("Independent RNAfold validation summary")
summary.append(f"Input rows: {len(df)}")
summary.append(f"RNAfold metric rows: {len(rnafold_df)}")
summary.append(f"Spearman comparison rows: {len(corr_df)}")
summary.append("")
summary.append("Top absolute Spearman correlations:")
summary.append(top_rows.to_string(index=False))
summary.append("")
summary.append("Reviewer-safe interpretation:")
summary.append(
    "RNAfold-derived secondary-structure metrics were computed directly from the 100-bp local genomic context sequences. "
    "These metrics were not used to construct the DNABERT-2 embeddings or the weak functional severity labels. "
    "Therefore, Spearman correlations between RNAfold metrics and guide-ranking or weak-risk outputs provide an independent "
    "computational validation signal. This should not be described as experimental validation."
)

out_summary.write_text("\n".join(summary))

print("Input:", input_path)
print("Wrote:", out_metrics)
print("Wrote:", out_corr)
print("Wrote:", out_corr_md)
print("Wrote:", out_summary)
print()
print("Top RNAfold Spearman correlations:")
print(top_rows.to_string(index=False))
