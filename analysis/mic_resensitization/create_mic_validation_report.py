from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[2]

candidate_file = root / "results_mic_resensitization/tables/mic_planning_candidate_guides.csv"
fold_summary_file = root / "results_mic_resensitization/tables/mic_fold_change_summary_template.csv"

out_report = root / "results_mic_resensitization/reports/mic_validation_planning_report.txt"
out_table = root / "results_mic_resensitization/tables/mic_validation_status_table.csv"

candidates = pd.read_csv(candidate_file)
fold_summary = pd.read_csv(fold_summary_file)

status_rows = []

for _, row in candidates.iterrows():
    gene = row["target_gene"]
    guide = row["guide_id"]
    role = row["mic_planning_role"]
    decision = row["final_decision"]

    if role == "Primary_MIC_candidate":
        validation_status = "Ready_for_future_MIC_planning"
    elif "Backup" in role or "backup" in role:
        validation_status = "Backup_for_future_MIC_planning"
    else:
        validation_status = "Caution_only_not_primary"

    status_rows.append(
        {
            "target_gene": gene,
            "guide_id": guide,
            "mic_planning_role": role,
            "requirement4_decision": decision,
            "experimental_MIC_status": "Not_done",
            "validation_status": validation_status,
        }
    )

status_table = pd.DataFrame(status_rows)
status_table.to_csv(out_table, index=False)

primary = status_table[status_table["mic_planning_role"] == "Primary_MIC_candidate"]
caution = status_table[status_table["validation_status"] == "Caution_only_not_primary"]

report = []
report.append("MIC Resensitization Validation Planning Report")
report.append("=" * 55)
report.append("")
report.append("Current status")
report.append("- MIC resensitization has not been experimentally performed.")
report.append("- No wet-lab MIC reduction, growth-curve rescue, plasmid curing, or cleavage validation result is claimed.")
report.append("- The current output is a computational planning module for future validation.")
report.append("")
report.append("Why this step was created")
report.append("- Requirements 1 to 4 completed the computational off-target safety framework.")
report.append("- MIC validation is the next biological validation step, but it remains pending.")
report.append("- This module prepares candidate guides, control conditions, and fold-change analysis templates.")
report.append("")
report.append("Primary computational candidates for future MIC testing")
for _, row in primary.iterrows():
    report.append(f"- {row['target_gene']}: {row['guide_id']} ({row['mic_planning_role']})")
report.append("")
report.append("Caution-only candidates")
for _, row in caution.iterrows():
    report.append(f"- {row['target_gene']}: {row['guide_id']} ({row['mic_planning_role']})")
report.append("")
report.append("Future MIC calculation plan")
report.append("- Add baseline MIC values into data_mic_resensitization/templates/mic_fold_change_input_template.csv")
report.append("- Add post-CRISPR MIC values into the same file.")
report.append("- Re-run: python analysis/mic_resensitization/calculate_mic_fold_change.py")
report.append("- The script will calculate fold_reduction and log2_fold_reduction.")
report.append("")
report.append("Interpretation rule for future data")
report.append("- >=4-fold MIC reduction: strong resensitization signal.")
report.append("- >=2-fold MIC reduction: moderate resensitization signal.")
report.append("- <2-fold MIC reduction: no clear resensitization signal.")
report.append("")
report.append("Manuscript-safe wording")
report.append(
    "MIC-resensitization experiments were not performed in the current study. "
    "Instead, the final FOTR-CRISPR guide set was translated into a validation-ready MIC planning framework, "
    "including prioritized candidate guides, proposed control conditions, antibiotic contexts, and a fold-change analysis template "
    "for future experimental measurements."
)
report.append("")
report.append("Output files")
report.append(f"- {out_table}")
report.append(f"- {out_report}")

out_report.write_text("\n".join(report))

print("MIC validation planning report created.")
print("Wrote:", out_table)
print("Wrote:", out_report)
