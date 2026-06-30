from pathlib import Path
import pandas as pd
import numpy as np
import re

ROOT = Path(__file__).resolve().parents[2]

OUT_TABLES = ROOT / "results_journal_extension/tables"
OUT_REPORTS = ROOT / "results_journal_extension/reports"
OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_REPORTS.mkdir(parents=True, exist_ok=True)

REQ4_FILE = ROOT / "results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/final_requirement4_complete/requirement4_final_integrated_conserved_offtarget_table_ALL20.csv"
REQ2_HITS_FILE = ROOT / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2_database_refined/requirement2_database_refined_hit_annotations.csv"
MIC_STATUS_FILE = ROOT / "results_mic_resensitization/tables/mic_validation_status_table.csv"
MIC_CANDIDATE_FILE = ROOT / "results_mic_resensitization/tables/mic_planning_candidate_guides.csv"

if not REQ4_FILE.exists():
    raise FileNotFoundError(f"Missing Requirement 4 final table: {REQ4_FILE}")

req4 = pd.read_csv(REQ4_FILE)

def get_col(df, preferred):
    for c in preferred:
        if c in df.columns:
            return c
    lower = {c.lower(): c for c in df.columns}
    for c in preferred:
        if c.lower() in lower:
            return lower[c.lower()]
    for col in df.columns:
        for c in preferred:
            if c.lower() in col.lower():
                return col
    return None

gene_col = get_col(req4, ["target_gene", "gene"])
guide_col = get_col(req4, ["guide_id"])
spacer_col = get_col(req4, ["spacer"])
hits_col = get_col(req4, ["total_offtarget_hit_rows", "total_offtarget_hits", "offtarget_hits"])
acc_col = get_col(req4, ["accessions_with_any_hit", "accessions_with_hits", "genomes_with_hits"])
frac_col = get_col(req4, ["fraction_accessions_with_any_hit"])
func_burden_col = get_col(req4, ["total_db_refined_burden", "db_refined_total_burden"])
mobile_burden_col = get_col(req4, ["requirement3_mobile_context_burden"])
cons_burden_col = get_col(req4, ["conserved_functional_burden_score"])
score_col = get_col(req4, ["final_req2_req3_req4_integrated_score", "requirement4_conservation_integrated_score"])
decision_col = get_col(req4, ["final_req4_integrated_decision", "requirement4_final_conservation_decision"])
rank_col = get_col(req4, ["final_req4_rank_within_gene"])

needed = {
    "gene_col": gene_col,
    "guide_col": guide_col,
    "spacer_col": spacer_col,
    "hits_col": hits_col,
    "acc_col": acc_col,
    "func_burden_col": func_burden_col,
    "mobile_burden_col": mobile_burden_col,
    "cons_burden_col": cons_burden_col,
    "decision_col": decision_col,
    "rank_col": rank_col,
}

missing = [k for k, v in needed.items() if v is None]
if missing:
    raise ValueError(f"Missing required columns: {missing}\nAvailable columns: {list(req4.columns)}")

# -------------------------
# Table 1: all-20 guide table
# -------------------------
all20_cols = [
    gene_col,
    guide_col,
    spacer_col,
    hits_col,
    acc_col,
    func_burden_col,
    mobile_burden_col,
    cons_burden_col,
    decision_col,
    rank_col,
]
if score_col:
    all20_cols.insert(-1, score_col)
if frac_col:
    all20_cols.insert(5, frac_col)

all20 = req4[all20_cols].copy()

rename_map = {
    gene_col: "target_gene",
    guide_col: "guide_id",
    spacer_col: "spacer",
    hits_col: "total_offtarget_hits",
    acc_col: "accessions_with_any_hit",
    func_burden_col: "functional_burden_score",
    mobile_burden_col: "mobile_context_burden_score",
    cons_burden_col: "conserved_offtarget_burden_score",
    decision_col: "final_integrated_decision",
    rank_col: "final_rank_within_gene",
}
if score_col:
    rename_map[score_col] = "final_integrated_score"
if frac_col:
    rename_map[frac_col] = "fraction_accessions_with_any_hit"

all20 = all20.rename(columns=rename_map)
all20 = all20.sort_values(["target_gene", "final_rank_within_gene"])

all20.to_csv(OUT_TABLES / "journal_table_all20_integrated_fotr_guides.csv", index=False)

# -------------------------
# Table 2: top two per gene
# -------------------------
top2 = all20[all20["final_rank_within_gene"] <= 2].copy()
top2.to_csv(OUT_TABLES / "journal_table_top2_guides_per_gene.csv", index=False)

# -------------------------
# Table 3: deprioritized guides
# -------------------------
deprioritized = all20[
    all20["final_integrated_decision"].astype(str).str.contains("DEPRIORITIZE|Deprioritize|high", case=False, na=False)
].copy()
deprioritized.to_csv(OUT_TABLES / "journal_table_deprioritized_guides.csv", index=False)

# -------------------------
# Table 4: requirement status table
# -------------------------
status_rows = [
    {
        "requirement": "Requirement 1",
        "analysis_layer": "Expanded Cas-OFFinder genome-wide off-target hit table",
        "journal_status": "Complete computationally",
        "main_evidence": "36-accession expanded panel; 1,997 parsed off-target hits; 20/20 guides audited",
        "manuscript_claim": "Real candidate off-target sites were generated and mapped to selected guides."
    },
    {
        "requirement": "Requirement 2",
        "analysis_layer": "Functional annotation of off-target hits",
        "journal_status": "Complete computationally",
        "main_evidence": "1,997/1,997 hits annotated using matched GFF features and database-refined dictionaries",
        "manuscript_claim": "Candidate off-targets were assigned genomic and functional context."
    },
    {
        "requirement": "Requirement 3",
        "analysis_layer": "Mobile-context, plasmid, SCCmec, integron, IS, transposon mapping",
        "journal_status": "Complete computationally",
        "main_evidence": "+/-10 kb neighborhood mapping plus direct evidence layers",
        "manuscript_claim": "Mobile-genetic context was integrated into off-target prioritization."
    },
    {
        "requirement": "Requirement 4",
        "analysis_layer": "Off-target conservation across strain genomes",
        "journal_status": "Complete computationally",
        "main_evidence": "Guide-level and site-level recurrence across 36 accessions",
        "manuscript_claim": "Recurring off-target burden was quantified across the expanded genome panel."
    },
    {
        "requirement": "MIC planning",
        "analysis_layer": "MIC-resensitization validation planning",
        "journal_status": "Planning complete; wet lab not performed",
        "main_evidence": "Candidate/control/fold-change templates created",
        "manuscript_claim": "The final guide set was translated into a validation-ready MIC planning framework."
    },
    {
        "requirement": "Wet-lab MIC validation",
        "analysis_layer": "Experimental MIC testing",
        "journal_status": "Not done",
        "main_evidence": "No experimental MIC values, growth curves, cleavage assays, or plasmid-curing data available",
        "manuscript_claim": "No wet-lab resensitization claim is made."
    },
]

status = pd.DataFrame(status_rows)
status.to_csv(OUT_TABLES / "journal_table_requirement_completion_status.csv", index=False)

# -------------------------
# Table 5: MIC planning table if available
# -------------------------
if MIC_STATUS_FILE.exists():
    mic_status = pd.read_csv(MIC_STATUS_FILE)
    mic_status.to_csv(OUT_TABLES / "journal_table_mic_validation_planning_status.csv", index=False)

if MIC_CANDIDATE_FILE.exists():
    mic_candidates = pd.read_csv(MIC_CANDIDATE_FILE)
    mic_candidates.to_csv(OUT_TABLES / "journal_table_mic_candidate_guides.csv", index=False)

# -------------------------
# Table 6: database-refined annotation validation summary
# -------------------------
db_summary_rows = []

if REQ2_HITS_FILE.exists():
    hits = pd.read_csv(REQ2_HITS_FILE)

    def count_signal(columns_keywords):
        cols = []
        for col in hits.columns:
            cl = col.lower()
            if any(k in cl for k in columns_keywords):
                cols.append(col)
        if not cols:
            return None, 0

        signal = pd.Series(False, index=hits.index)
        for col in cols:
            s = hits[col]
            if s.dtype == bool:
                signal = signal | s.fillna(False)
            else:
                sval = s.astype(str).str.lower()
                signal = signal | sval.isin(["true", "1", "yes", "y"])
                signal = signal | sval.str.contains("amr|resistance|virulence|essential|mobile|plasmid|integron|transposon|insertion", na=False)
        return cols, int(signal.sum())

    categories = [
        ("AMR_or_resistance_signal", ["amr", "resistance", "antibiotic"]),
        ("Virulence_signal", ["virulence", "vfdb"]),
        ("Essential_like_signal", ["essential", "deg"]),
        ("Mobile_or_plasmid_signal", ["mobile", "plasmid", "integron", "transposon", "insertion", "is_"]),
        ("High_confidence_functional_signal", ["high_confidence"]),
    ]

    for label, keys in categories:
        cols, n = count_signal(keys)
        db_summary_rows.append({
            "database_refined_signal": label,
            "supporting_columns_used": "; ".join(cols) if cols else "No direct column found",
            "hit_rows_with_signal": n,
            "total_hit_rows": len(hits),
            "fraction_of_hits": round(n / len(hits), 4) if len(hits) else 0
        })

    db_summary = pd.DataFrame(db_summary_rows)
    db_summary.to_csv(OUT_TABLES / "journal_table_database_refined_annotation_validation_summary.csv", index=False)

# -------------------------
# Journal methods/results text
# -------------------------
report = []
report.append("Journal Extension Computational Analysis Summary")
report.append("=" * 60)
report.append("")
report.append("Scope:")
report.append("This journal extension is fully computational. Wet-lab MIC-resensitization validation was not performed.")
report.append("")
report.append("Journal-ready tables created:")
for f in sorted(OUT_TABLES.glob("journal_table_*.csv")):
    report.append(f"- {f}")
report.append("")
report.append("Primary computational conclusion:")
report.append("The final integrated FOTR-CRISPR ranking identifies blaKPC_riskaware_top5 and blaNDM1_riskaware_top1 as the cleanest computational candidates under the current 36-accession off-target, functional, mobile-context, and conservation framework.")
report.append("")
report.append("Safe MIC wording:")
report.append("MIC-resensitization experiments were not performed. The current work provides a validation-ready MIC planning framework only.")
report.append("")
report.append("Suggested journal claim:")
report.append("This study presents a computational FOTR-CRISPR framework that integrates genome-wide off-target discovery, functional annotation, mobile-genetic context, and strain-level off-target conservation to prioritize CRISPR-Cas9 guide RNAs against AMR genes.")

(OUT_REPORTS / "journal_extension_table_generation_summary.txt").write_text("\n".join(report))

print("Journal-ready tables created.")
for f in sorted(OUT_TABLES.glob("journal_table_*.csv")):
    print("Wrote:", f)
print("Wrote:", OUT_REPORTS / "journal_extension_table_generation_summary.txt")
