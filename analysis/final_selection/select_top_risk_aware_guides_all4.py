import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results_mlcb" / "tables" / "functional_offtarget_risk_ranked_guides.csv"

OUTPUT_DIR = PROJECT_ROOT / "results_final_selection"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_N_PER_GENE = 5

CORE_GENES = ["blaKPC", "blaNDM1", "mcr1", "mecA"]


def main():
    print("Selecting top risk-aware guides for all 4 core AMR genes...")
    print(f"Input file: {INPUT_FILE}")

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Missing input file: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded ranked guide table: {df.shape}")

    required_cols = [
        "mlcb_gene",
        "mlcb_spacer",
        "mlcb_base_score",
        "annotation_derived_functional_risk_score",
        "predicted_high_functional_offtarget_risk_probability",
        "functional_risk_penalized_score",
        "risk_aware_rank",
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Keep only the 4 original AMR genes
    df = df[df["mlcb_gene"].isin(CORE_GENES)].copy()
    print(f"Rows after filtering to core genes: {df.shape}")

    # Remove duplicate gene-spacer rows.
    # Your current output has duplicates because the same spacer may appear with multiple PAM/context records.
    df = df.sort_values(
        by=[
            "mlcb_gene",
            "functional_risk_penalized_score",
            "predicted_high_functional_offtarget_risk_probability",
        ],
        ascending=[True, False, True],
    )

    dedup = df.drop_duplicates(subset=["mlcb_gene", "mlcb_spacer"], keep="first").copy()
    print(f"Rows after gene-spacer deduplication: {dedup.shape}")

    # Re-rank within each gene
    dedup["final_gene_rank"] = (
        dedup.groupby("mlcb_gene")["functional_risk_penalized_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )

    # Keep top N per gene
    top = dedup[dedup["final_gene_rank"] <= TOP_N_PER_GENE].copy()
    top = top.sort_values(["mlcb_gene", "final_gene_rank"])

    # Add guide IDs for downstream Cas-OFFinder / MIC planning
    top["final_guide_id"] = (
        top["mlcb_gene"].astype(str)
        + "_riskaware_top"
        + top["final_gene_rank"].astype(str)
    )

    selected_cols = [
        "final_guide_id",
        "mlcb_gene",
        "mlcb_spacer",
        "final_gene_rank",
        "mlcb_base_score",
        "annotation_derived_functional_risk_score",
        "predicted_high_functional_offtarget_risk_probability",
        "functional_risk_penalized_score",
    ]

    top_out = top[selected_cols].copy()

    all_outfile = OUTPUT_DIR / "all_core_guides_deduplicated_risk_aware_ranked.csv"
    top_outfile = OUTPUT_DIR / "top5_per_gene_risk_aware_guides.csv"
    cas_input_file = OUTPUT_DIR / "cas_offinder_input_guides_top5_per_gene.tsv"
    fasta_file = OUTPUT_DIR / "top5_per_gene_risk_aware_guides.fasta"
    summary_file = OUTPUT_DIR / "top5_per_gene_risk_aware_summary.txt"

    dedup.to_csv(all_outfile, index=False)
    top_out.to_csv(top_outfile, index=False)

    # TSV for Cas-OFFinder preparation
    cas_df = top_out[["final_guide_id", "mlcb_gene", "mlcb_spacer"]].copy()
    cas_df.to_csv(cas_input_file, sep="\t", index=False)

    # FASTA version
    with open(fasta_file, "w") as f:
        for _, row in top_out.iterrows():
            f.write(f">{row['final_guide_id']}|{row['mlcb_gene']}|rank_{row['final_gene_rank']}\n")
            f.write(f"{row['mlcb_spacer']}\n")

    with open(summary_file, "w") as f:
        f.write("Top Risk-Aware Guide Selection Summary\n")
        f.write("=====================================\n\n")
        f.write(f"Input ranked table: {INPUT_FILE}\n")
        f.write(f"Original filtered rows for 4 genes: {df.shape[0]}\n")
        f.write(f"Deduplicated gene-spacer rows: {dedup.shape[0]}\n")
        f.write(f"Top guides selected per gene: {TOP_N_PER_GENE}\n\n")

        for gene in CORE_GENES:
            gene_top = top_out[top_out["mlcb_gene"] == gene]
            f.write(f"{gene}: {len(gene_top)} selected guides\n")
            for _, row in gene_top.iterrows():
                f.write(
                    f"  Rank {row['final_gene_rank']}: {row['mlcb_spacer']} "
                    f"| penalized_score={row['functional_risk_penalized_score']:.6f} "
                    f"| predicted_risk={row['predicted_high_functional_offtarget_risk_probability']:.6f}\n"
                )
            f.write("\n")

    print("\nWrote output files:")
    print(f"- {all_outfile}")
    print(f"- {top_outfile}")
    print(f"- {cas_input_file}")
    print(f"- {fasta_file}")
    print(f"- {summary_file}")

    print("\nTop selected guides:")
    print(top_out.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()