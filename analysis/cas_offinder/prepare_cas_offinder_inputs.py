import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

GUIDE_FILE = PROJECT_ROOT / "results_final_selection" / "top5_per_gene_risk_aware_guides.csv"

INPUT_DIR = PROJECT_ROOT / "data_cas_offinder" / "input"
INPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = INPUT_DIR / "top5_per_gene_cas_offinder_queries.tsv"
SUMMARY_FILE = INPUT_DIR / "cas_offinder_input_summary.txt"

CORE_GENES = ["blaKPC", "blaNDM1", "mcr1", "mecA"]


def main():
    print("Preparing Cas-OFFinder query table for top risk-aware guides...")

    if not GUIDE_FILE.exists():
        raise FileNotFoundError(f"Missing guide file: {GUIDE_FILE}")

    df = pd.read_csv(GUIDE_FILE)
    print(f"Loaded guide table: {df.shape}")

    required = ["final_guide_id", "mlcb_gene", "mlcb_spacer", "final_gene_rank"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df[df["mlcb_gene"].isin(CORE_GENES)].copy()

    # Cas-OFFinder commonly searches target pattern: 20nt guide + NGG PAM.
    # We do not know the exact genomic PAM here, so use NNN-compatible PAM pattern separately.
    # For SpCas9, the guide sequence is 20 nt and PAM is NGG.
    df["cas9_pattern"] = df["mlcb_spacer"] + "NGG"
    df["allowed_mismatches"] = 4

    out = df[
        [
            "final_guide_id",
            "mlcb_gene",
            "final_gene_rank",
            "mlcb_spacer",
            "cas9_pattern",
            "allowed_mismatches",
            "functional_risk_penalized_score",
            "predicted_high_functional_offtarget_risk_probability",
        ]
    ].copy()

    out.to_csv(OUTPUT_FILE, sep="\t", index=False)

    with open(SUMMARY_FILE, "w") as f:
        f.write("Cas-OFFinder Input Preparation Summary\n")
        f.write("=====================================\n\n")
        f.write(f"Input guide file: {GUIDE_FILE}\n")
        f.write(f"Total guides prepared: {out.shape[0]}\n\n")
        for gene in CORE_GENES:
            sub = out[out["mlcb_gene"] == gene]
            f.write(f"{gene}: {sub.shape[0]} guides\n")
            for _, row in sub.iterrows():
                f.write(
                    f"  {row['final_guide_id']}: {row['mlcb_spacer']} + NGG "
                    f"| rank={row['final_gene_rank']} "
                    f"| mismatches={row['allowed_mismatches']}\n"
                )
            f.write("\n")

    print("\nWrote:")
    print(f"- {OUTPUT_FILE}")
    print(f"- {SUMMARY_FILE}")

    print("\nPreview:")
    print(out.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()