#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(".").resolve()

GUIDE_FILE = PROJECT_ROOT / "results_final_selection/cas_offinder_input_guides_top5_per_gene.tsv"
HIT_FILE = PROJECT_ROOT / "results_cas_offinder/expanded_panel/parsed/expanded_casoffinder_all_genome_offtarget_hits_parsed.csv"
MANIFEST_FILE = PROJECT_ROOT / "data_cas_offinder/expanded_panel/metadata/expanded_casoffinder_run_manifest.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/qc"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_AUDIT = OUT_DIR / "expanded_casoffinder_guide_completeness_audit.csv"
OUT_ZERO = OUT_DIR / "expanded_casoffinder_zero_hit_guides.csv"
OUT_SUMMARY = OUT_DIR / "expanded_casoffinder_qc_summary.txt"


def load_guides():
    guides = pd.read_csv(GUIDE_FILE, sep="\t")
    guides.columns = [c.strip() for c in guides.columns]

    spacer_col = "mlcb_spacer" if "mlcb_spacer" in guides.columns else "spacer"
    guide_col = "final_guide_id" if "final_guide_id" in guides.columns else "guide_id"
    gene_col = "mlcb_gene" if "mlcb_gene" in guides.columns else "gene"

    guides = guides.rename(columns={
        spacer_col: "spacer",
        guide_col: "guide_id",
        gene_col: "target_gene"
    })

    guides["spacer"] = guides["spacer"].astype(str).str.upper().str.strip()
    guides["guide_id"] = guides["guide_id"].astype(str)
    guides["target_gene"] = guides["target_gene"].astype(str)

    return guides[["guide_id", "target_gene", "spacer"]].drop_duplicates()


def main():
    guides = load_guides()
    hits = pd.read_csv(HIT_FILE)
    manifest = pd.read_csv(MANIFEST_FILE)

    guide_summary = (
        hits.groupby(["guide_id", "target_gene", "spacer"], dropna=False)
        .agg(
            total_expanded_hits=("guide_id", "size"),
            genomes_with_hits=("accession", "nunique"),
            min_mismatch=("mismatch_count", "min"),
            max_mismatch=("mismatch_count", "max"),
            mismatch_1_hits=("mismatch_count", lambda x: (x == 1).sum()),
            mismatch_2_hits=("mismatch_count", lambda x: (x == 2).sum()),
            mismatch_3_hits=("mismatch_count", lambda x: (x == 3).sum()),
            mismatch_4_hits=("mismatch_count", lambda x: (x == 4).sum()),
        )
        .reset_index()
    )

    audit = guides.merge(
        guide_summary,
        on=["guide_id", "target_gene", "spacer"],
        how="left"
    )

    fill_cols = [
        "total_expanded_hits",
        "genomes_with_hits",
        "mismatch_1_hits",
        "mismatch_2_hits",
        "mismatch_3_hits",
        "mismatch_4_hits",
    ]

    for col in fill_cols:
        audit[col] = audit[col].fillna(0).astype(int)

    audit["expanded_hit_status"] = audit["total_expanded_hits"].apply(
        lambda x: "HIT_POSITIVE" if x > 0 else "ZERO_HIT_IN_EXPANDED_PANEL"
    )

    zero = audit[audit["expanded_hit_status"] == "ZERO_HIT_IN_EXPANDED_PANEL"].copy()

    audit.to_csv(OUT_AUDIT, index=False)
    zero.to_csv(OUT_ZERO, index=False)

    expected_guides = guides["guide_id"].nunique()
    observed_guides = audit["guide_id"].nunique()
    hit_positive_guides = (audit["total_expanded_hits"] > 0).sum()
    zero_hit_guides = (audit["total_expanded_hits"] == 0).sum()

    invalid_mismatches = hits[~hits["mismatch_count"].isin([0, 1, 2, 3, 4])].shape[0]
    invalid_positions = hits[hits["position_1based"] <= 0].shape[0]
    invalid_strands = hits[~hits["strand"].isin(["+", "-"])].shape[0]

    summary = []
    summary.append("Expanded Cas-OFFinder Completeness QC")
    summary.append("====================================")
    summary.append(f"Expected selected guides: {expected_guides}")
    summary.append(f"Observed selected guides in audit: {observed_guides}")
    summary.append(f"Expanded genome configs/runs expected: {manifest.shape[0]}")
    summary.append(f"Unique accessions in hit table: {hits['accession'].nunique()}")
    summary.append(f"Total expanded off-target hits: {hits.shape[0]}")
    summary.append(f"Hit-positive guides: {hit_positive_guides}")
    summary.append(f"Zero-hit guides in expanded panel: {zero_hit_guides}")
    summary.append(f"Invalid mismatch rows: {invalid_mismatches}")
    summary.append(f"Invalid position rows: {invalid_positions}")
    summary.append(f"Invalid strand rows: {invalid_strands}")
    summary.append("")
    summary.append("Hit count by target gene:")
    summary.append(str(audit.groupby("target_gene")["total_expanded_hits"].sum().sort_values(ascending=False)))
    summary.append("")
    summary.append("Top guides by expanded hit count:")
    summary.append(str(audit.sort_values("total_expanded_hits", ascending=False)[["guide_id", "target_gene", "total_expanded_hits", "genomes_with_hits"]].head(10).to_string(index=False)))
    summary.append("")
    summary.append("Zero-hit guides:")
    if zero.empty:
        summary.append("None")
    else:
        summary.append(str(zero[["guide_id", "target_gene", "spacer"]].to_string(index=False)))
    summary.append("")
    summary.append(f"Audit table: {OUT_AUDIT}")
    summary.append(f"Zero-hit table: {OUT_ZERO}")

    qc_pass = (
        expected_guides == 20
        and observed_guides == 20
        and manifest.shape[0] == 36
        and hits["accession"].nunique() <= 36
        and invalid_mismatches == 0
        and invalid_positions == 0
        and invalid_strands == 0
    )

    summary.append("QC status: " + ("PASS" if qc_pass else "CHECK"))

    OUT_SUMMARY.write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
