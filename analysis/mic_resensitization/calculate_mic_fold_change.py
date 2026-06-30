import pandas as pd
import numpy as np
from pathlib import Path

root = Path(__file__).resolve().parents[2]

input_file = root / "data_mic_resensitization/templates/mic_fold_change_input_template.csv"
out_dir = root / "results_mic_resensitization/tables"
report_dir = root / "results_mic_resensitization/reports"

out_dir.mkdir(parents=True, exist_ok=True)
report_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

# Convert MIC columns to numeric if future values are added.
df["baseline_MIC_no_CRISPR"] = pd.to_numeric(df["baseline_MIC_no_CRISPR"], errors="coerce")
df["post_CRISPR_MIC"] = pd.to_numeric(df["post_CRISPR_MIC"], errors="coerce")

# Calculate only where both values exist.
df["fold_reduction"] = df["baseline_MIC_no_CRISPR"] / df["post_CRISPR_MIC"]
df["log2_fold_reduction"] = np.log2(df["fold_reduction"])

def call_resensitization(row):
    if pd.isna(row["fold_reduction"]):
        return row["resensitization_call"]
    if row["fold_reduction"] >= 4:
        return "Strong_resensitization_signal"
    if row["fold_reduction"] >= 2:
        return "Moderate_resensitization_signal"
    return "No_clear_resensitization_signal"

df["resensitization_call"] = df.apply(call_resensitization, axis=1)

df.to_csv(out_dir / "mic_fold_change_results_template.csv", index=False)

summary = df.groupby(["target_gene", "guide_id", "antibiotic"], dropna=False).agg(
    replicates=("replicate", "count"),
    mean_baseline_MIC=("baseline_MIC_no_CRISPR", "mean"),
    mean_post_CRISPR_MIC=("post_CRISPR_MIC", "mean"),
    mean_fold_reduction=("fold_reduction", "mean"),
    mean_log2_fold_reduction=("log2_fold_reduction", "mean"),
).reset_index()

summary.to_csv(out_dir / "mic_fold_change_summary_template.csv", index=False)

report = """MIC Fold-Change Analysis Template

Status:
No experimental MIC values have been added yet.

Purpose:
This script is ready to calculate MIC fold-reduction after future wet-lab MIC measurements are available.

Calculation:
fold_reduction = baseline_MIC_no_CRISPR / post_CRISPR_MIC
log2_fold_reduction = log2(fold_reduction)

Interpretation rule:
>=4-fold reduction: strong resensitization signal
>=2-fold reduction: moderate resensitization signal
<2-fold reduction: no clear resensitization signal

Important:
Blank MIC values mean experimental validation is still pending.
"""

(report_dir / "mic_fold_change_analysis_template_summary.txt").write_text(report)

print("MIC fold-change template analysis completed.")
print("Wrote:", out_dir / "mic_fold_change_results_template.csv")
print("Wrote:", out_dir / "mic_fold_change_summary_template.csv")
print("Wrote:", report_dir / "mic_fold_change_analysis_template_summary.txt")
