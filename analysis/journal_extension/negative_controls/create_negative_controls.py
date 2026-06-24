from pathlib import Path
import random
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[3]

INPUT_FOTR_V2 = ROOT / "results_fotr_v2" / "all_guides_fotr_v2_functional_context_recommended.csv"

OUT_DIR = ROOT / "results_journal_extension" / "negative_controls"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "negative_control_guides.csv"
OUT_SUMMARY = OUT_DIR / "negative_control_summary.txt"

random.seed(42)

DNA = ["A", "C", "G", "T"]


def random_spacer(length=20):
    return "".join(random.choice(DNA) for _ in range(length))


def gc_content(seq):
    seq = str(seq).upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def make_controls(real_df, n_per_gene=25):
    rows = []

    genes = sorted(real_df["gene"].dropna().unique())

    for gene in genes:
        real_gene = real_df[real_df["gene"] == gene].copy()

        # Control 1: fully random 20-nt guides
        for i in range(n_per_gene):
            spacer = random_spacer()
            rows.append({
                "gene": gene,
                "control_type": "random_20nt_sequence",
                "spacer": spacer,
                "pam": "NGG",
                "gc_content": gc_content(spacer),
                "control_reason": "Random sequence not selected from AMR target gene.",
                "expected_priority": "Low"
            })

        # Control 2: GC-extreme guides
        for i in range(max(5, n_per_gene // 5)):
            spacer_high_gc = "".join(random.choice(["G", "C"]) for _ in range(20))
            spacer_low_gc = "".join(random.choice(["A", "T"]) for _ in range(20))

            rows.append({
                "gene": gene,
                "control_type": "extreme_high_gc_sequence",
                "spacer": spacer_high_gc,
                "pam": "NGG",
                "gc_content": gc_content(spacer_high_gc),
                "control_reason": "Extreme GC content can reduce practical guide suitability.",
                "expected_priority": "Low"
            })

            rows.append({
                "gene": gene,
                "control_type": "extreme_low_gc_sequence",
                "spacer": spacer_low_gc,
                "pam": "NGG",
                "gc_content": gc_content(spacer_low_gc),
                "control_reason": "Extreme AT content can reduce practical guide suitability.",
                "expected_priority": "Low"
            })

        # Control 3: weakest real candidates by FOTR v2
        if "fotr_v2_priority_score" in real_gene.columns:
            weakest = real_gene.sort_values("fotr_v2_priority_score", ascending=True).head(5)
            for _, r in weakest.iterrows():
                rows.append({
                    "gene": gene,
                    "control_type": "lowest_fotr_v2_real_candidate",
                    "spacer": r.get("spacer", ""),
                    "pam": r.get("pam", ""),
                    "gc_content": r.get("gc_content", np.nan),
                    "control_reason": "Real candidate with lowest FOTR v2 priority score for this gene.",
                    "expected_priority": "Low",
                    "fotr_v2_priority_score": r.get("fotr_v2_priority_score", np.nan),
                    "rna_accessibility_class": r.get("rna_accessibility_class", ""),
                    "fotr_v2_recommendation_status": r.get("fotr_v2_recommendation_status", "")
                })

        # Control 4: structurally risky real candidates
        if "rna_accessibility_class" in real_gene.columns:
            risky = real_gene[real_gene["rna_accessibility_class"].astype(str).str.contains("Structurally_Risky", na=False)].head(5)
            for _, r in risky.iterrows():
                rows.append({
                    "gene": gene,
                    "control_type": "structurally_risky_real_candidate",
                    "spacer": r.get("spacer", ""),
                    "pam": r.get("pam", ""),
                    "gc_content": r.get("gc_content", np.nan),
                    "control_reason": "Real candidate classified as structurally risky by RNA accessibility screen.",
                    "expected_priority": "Low_or_filtered",
                    "fotr_v2_priority_score": r.get("fotr_v2_priority_score", np.nan),
                    "rna_accessibility_class": r.get("rna_accessibility_class", ""),
                    "fotr_v2_recommendation_status": r.get("fotr_v2_recommendation_status", "")
                })

    return pd.DataFrame(rows)


def main():
    print("Reading FOTR v2 recommended table...")
    df = pd.read_csv(INPUT_FOTR_V2)

    controls = make_controls(df, n_per_gene=25)

    # Add simple control suitability flags
    controls["gc_flag"] = np.where(
        (controls["gc_content"] < 30) | (controls["gc_content"] > 70),
        "GC_extreme",
        "GC_reasonable"
    )

    controls["negative_control_label"] = "Negative_control"

    controls.to_csv(OUT_TABLE, index=False)

    summary = []
    summary.append("Negative Control Summary")
    summary.append("========================")
    summary.append("")
    summary.append(f"Input FOTR v2 table: {INPUT_FOTR_V2}")
    summary.append(f"Output table: {OUT_TABLE}")
    summary.append("")
    summary.append(f"Total negative-control rows: {len(controls)}")
    summary.append(f"Genes represented: {', '.join(sorted(controls['gene'].dropna().astype(str).unique()))}")
    summary.append("")
    summary.append("Control type counts:")
    for k, v in controls["control_type"].value_counts().items():
        summary.append(f"- {k}: {v}")
    summary.append("")
    summary.append("GC flag counts:")
    for k, v in controls["gc_flag"].value_counts().items():
        summary.append(f"- {k}: {v}")
    summary.append("")
    summary.append("Interpretation:")
    summary.append(
        "These controls provide expected-low-priority examples, including random 20-nt sequences, "
        "GC-extreme synthetic sequences, lowest-ranked real candidates, and structurally risky real candidates. "
        "They are intended as computational sanity checks rather than experimental negative controls."
    )

    OUT_SUMMARY.write_text("\n".join(summary))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
