import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

candidate_file = PROJECT_ROOT / "data_mic_resensitization/templates/mic_candidate_guides_template.csv"
design_file = PROJECT_ROOT / "data_mic_resensitization/templates/mic_experimental_design_template.csv"

out_dir = PROJECT_ROOT / "results_mic_resensitization/tables"
report_dir = PROJECT_ROOT / "results_mic_resensitization/reports"

out_dir.mkdir(parents=True, exist_ok=True)
report_dir.mkdir(parents=True, exist_ok=True)

candidates = pd.read_csv(candidate_file)
design = pd.read_csv(design_file)

# Table 1: clean candidate summary
candidate_summary = candidates[
    [
        "target_gene",
        "guide_id",
        "spacer",
        "requirement4_rank",
        "offtarget_hits",
        "accessions_with_hits",
        "conservation_label",
        "final_decision",
        "mic_planning_role",
    ]
].copy()

candidate_summary.to_csv(
    out_dir / "mic_planning_candidate_guides.csv",
    index=False,
)

# Table 2: only primary/backup practical candidates
priority_candidates = candidate_summary[
    candidate_summary["mic_planning_role"].isin(
        ["Primary_MIC_candidate", "Backup_MIC_candidate", "Backup_or_caution_candidate"]
    )
].copy()

priority_candidates.to_csv(
    out_dir / "mic_priority_and_backup_candidates.csv",
    index=False,
)

# Table 3: MIC experimental design table
design.to_csv(
    out_dir / "mic_resensitization_experimental_design_template.csv",
    index=False,
)

# Table 4: manuscript-safe interpretation table
interpretation_rows = []

for _, row in candidate_summary.iterrows():
    gene = row["target_gene"]
    guide = row["guide_id"]
    role = row["mic_planning_role"]
    decision = row["final_decision"]
    hits = row["offtarget_hits"]
    acc = row["accessions_with_hits"]
    conservation = row["conservation_label"]

    if role == "Primary_MIC_candidate":
        interpretation = (
            f"{guide} is a primary computational candidate for future MIC-resensitization testing "
            f"because it showed {hits} off-target hits across {acc} accessions and was classified as {conservation}."
        )
    elif role == "Backup_MIC_candidate":
        interpretation = (
            f"{guide} is a backup computational candidate for future MIC-resensitization testing. "
            f"It has limited off-target recurrence but should be interpreted after the primary candidate."
        )
    elif "Caution" in role or "caution" in role:
        interpretation = (
            f"{guide} should be treated as a caution-only candidate because the current conservation-integrated "
            f"analysis indicates elevated off-target concern."
        )
    else:
        interpretation = (
            f"{guide} is retained for planning but requires careful review before experimental prioritization."
        )

    interpretation_rows.append(
        {
            "target_gene": gene,
            "guide_id": guide,
            "mic_planning_role": role,
            "requirement4_final_decision": decision,
            "manuscript_safe_interpretation": interpretation,
        }
    )

interpretation_table = pd.DataFrame(interpretation_rows)
interpretation_table.to_csv(
    out_dir / "mic_candidate_manuscript_safe_interpretation.csv",
    index=False,
)

# Text report
report_text = []
report_text.append("MIC Resensitization Planning Summary")
report_text.append("=" * 45)
report_text.append("")
report_text.append("Status: MIC resensitization has NOT been experimentally performed.")
report_text.append("Purpose: This module creates a computational planning table for future MIC validation.")
report_text.append("")
report_text.append("Input files:")
report_text.append(f"- {candidate_file}")
report_text.append(f"- {design_file}")
report_text.append("")
report_text.append("Output files:")
report_text.append(f"- {out_dir / 'mic_planning_candidate_guides.csv'}")
report_text.append(f"- {out_dir / 'mic_priority_and_backup_candidates.csv'}")
report_text.append(f"- {out_dir / 'mic_resensitization_experimental_design_template.csv'}")
report_text.append(f"- {out_dir / 'mic_candidate_manuscript_safe_interpretation.csv'}")
report_text.append("")
report_text.append("Key interpretation:")
report_text.append("- blaKPC_riskaware_top5 is the cleanest computational MIC-planning candidate.")
report_text.append("- blaNDM1_riskaware_top1 remains a strong low-recurrence computational candidate.")
report_text.append("- mcr1 and mecA candidates should be described cautiously because conserved off-target concern remains.")
report_text.append("")
report_text.append("Manuscript-safe wording:")
report_text.append(
    "MIC-resensitization has not yet been experimentally validated. "
    "The current analysis identifies computationally prioritized candidate guides and control conditions "
    "for future MIC testing."
)

(report_dir / "mic_resensitization_planning_summary.txt").write_text(
    "\n".join(report_text)
)

print("MIC planning tables created successfully.")
print(f"Wrote: {out_dir / 'mic_planning_candidate_guides.csv'}")
print(f"Wrote: {out_dir / 'mic_priority_and_backup_candidates.csv'}")
print(f"Wrote: {out_dir / 'mic_resensitization_experimental_design_template.csv'}")
print(f"Wrote: {out_dir / 'mic_candidate_manuscript_safe_interpretation.csv'}")
print(f"Wrote: {report_dir / 'mic_resensitization_planning_summary.txt'}")
