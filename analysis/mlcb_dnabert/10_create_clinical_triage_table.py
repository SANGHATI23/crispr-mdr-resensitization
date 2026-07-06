from pathlib import Path
import pandas as pd
import numpy as np

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
    "offtarget_penalty",
    "offtarget_hits_1mm",
    "offtarget_hits_2mm",
    "offtarget_hits_3mm",
    "number_of_reference_matches",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

work = df.copy()

numeric_cols = [
    "mlcb_base_score",
    "weak_functional_severity_label",
    "offtarget_penalty",
    "offtarget_hits_1mm",
    "offtarget_hits_2mm",
    "offtarget_hits_3mm",
    "number_of_reference_matches",
]

for c in numeric_cols:
    work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0)

def minmax(s):
    s = pd.to_numeric(s, errors="coerce").fillna(0)
    if s.max() == s.min():
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.min()) / (s.max() - s.min())

work["base_score_norm"] = minmax(work["mlcb_base_score"])

work["offtarget_hit_burden"] = (
    3.0 * work["offtarget_hits_1mm"]
    + 2.0 * work["offtarget_hits_2mm"]
    + 1.0 * work["offtarget_hits_3mm"]
)

work["risk_burden_score"] = (
    0.40 * minmax(work["offtarget_penalty"])
    + 0.25 * minmax(work["offtarget_hit_burden"])
    + 0.20 * minmax(work["weak_functional_severity_label"])
    + 0.15 * minmax(work["number_of_reference_matches"])
)

work["risk_aware_score"] = work["base_score_norm"] * (1.0 - work["risk_burden_score"])
work["estimated_high_risk_probability"] = minmax(work["risk_burden_score"])

rows = []

for gene, g in work.groupby("mlcb_gene"):
    g = g.copy()

    fotr_top = g.sort_values("mlcb_base_score", ascending=False).iloc[0]
    risk_top = g.sort_values("risk_aware_score", ascending=False).iloc[0]

    diverges = str(fotr_top["mlcb_spacer"]) != str(risk_top["mlcb_spacer"])

    if diverges:
        decision = "Switch from FOTR top-1 to risk-aware top-1 before wet-lab validation."
        action = (
            "Advance the risk-aware top guide first; keep the original FOTR top guide "
            "as backup; validate cleavage, resistance-gene disruption, resensitization, "
            "and off-target burden."
        )
    else:
        decision = "FOTR and risk-aware top-1 agree."
        action = (
            "Advance the same guide with higher confidence, while still validating "
            "cleavage, delivery, resensitization, and off-target burden experimentally."
        )

    rows.append({
        "AMR_gene": gene,

        "FOTR_top1_spacer": fotr_top["mlcb_spacer"],
        "FOTR_top1_base_score": round(float(fotr_top["mlcb_base_score"]), 4),
        "FOTR_top1_weak_severity_label": int(fotr_top["weak_functional_severity_label"]),
        "FOTR_top1_estimated_high_risk_probability": round(float(fotr_top["estimated_high_risk_probability"]), 4),
        "FOTR_top1_risk_aware_score": round(float(fotr_top["risk_aware_score"]), 4),

        "Risk_aware_top1_spacer": risk_top["mlcb_spacer"],
        "Risk_aware_top1_base_score": round(float(risk_top["mlcb_base_score"]), 4),
        "Risk_aware_top1_weak_severity_label": int(risk_top["weak_functional_severity_label"]),
        "Risk_aware_top1_estimated_high_risk_probability": round(float(risk_top["estimated_high_risk_probability"]), 4),
        "Risk_aware_top1_risk_aware_score": round(float(risk_top["risk_aware_score"]), 4),

        "Top1_diverges": diverges,
        "Clinical_triage_decision": decision,
        "Wet_lab_followup_action": action,
    })

triage = pd.DataFrame(rows)

triage["sort_key"] = triage["AMR_gene"].astype(str).str.lower().ne("mcr1").astype(int)
triage = triage.sort_values(["sort_key", "AMR_gene"]).drop(columns=["sort_key"])

out_csv = TABLE_DIR / "clinical_triage_decision_support_table.csv"
out_md = TABLE_DIR / "clinical_triage_decision_support_table.md"
out_txt = TABLE_DIR / "clinical_triage_decision_support_paragraph.txt"
ranked_out = TABLE_DIR / "clinical_triage_risk_aware_ranked_guides.csv"

triage.to_csv(out_csv, index=False)
triage.to_markdown(out_md, index=False)

work.sort_values(["mlcb_gene", "risk_aware_score"], ascending=[True, False]).to_csv(ranked_out, index=False)

paragraph = """
Clinical triage interpretation. The risk-aware ranking changes the practical wet-lab decision when the guide with the highest original FOTR score is not the same as the guide with the highest risk-aware score. In that setting, the workflow does not claim that one guide is experimentally safe; instead, it changes the order in which candidates should be validated. A wet-lab team would advance the risk-aware top-ranked guide first, retain the original FOTR top guide as a backup, and prioritize experimental checks for cleavage efficiency, resistance-gene disruption, bacterial resensitization, and off-target burden. This table therefore demonstrates the clinical utility of the computational layer: it converts model output into an auditable candidate-selection decision rather than only reporting classification metrics.
""".strip()

out_txt.write_text(paragraph)

print("Input:", input_path)
print("Wrote:", out_csv)
print("Wrote:", out_md)
print("Wrote:", out_txt)
print("Wrote:", ranked_out)
print()
print("Clinical triage table:")
print(triage.to_string(index=False))
print()
print("Divergence summary:")
print(triage[["AMR_gene", "Top1_diverges", "Clinical_triage_decision"]].to_string(index=False))
