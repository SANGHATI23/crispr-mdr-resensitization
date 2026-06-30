import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

GUIDE_FILE = PROJECT_ROOT / "results_final_selection" / "top5_per_gene_risk_aware_guides.csv"

GENOME_DIR = PROJECT_ROOT / "data_cas_offinder" / "genomes"
CONFIG_DIR = PROJECT_ROOT / "data_cas_offinder" / "configs"
CONFIG_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_FILE = CONFIG_DIR / "cas_offinder_config_summary.txt"

# SpCas9: 20-nt spacer + NGG PAM.
# Cas-OFFinder pattern length = 23.
CAS_OFFINDER_PATTERN = "NNNNNNNNNNNNNNNNNNNNNGG"

# We will search up to 4 mismatches.
ALLOWED_MISMATCHES = 4


def find_genome_fastas():
    fasta_extensions = ["*.fa", "*.fasta", "*.fna"]
    files = []
    for ext in fasta_extensions:
        files.extend(sorted(GENOME_DIR.glob(ext)))
    return files


def main():
    print("Creating Cas-OFFinder config files...")

    if not GUIDE_FILE.exists():
        raise FileNotFoundError(f"Missing guide file: {GUIDE_FILE}")

    guides = pd.read_csv(GUIDE_FILE)
    print(f"Loaded guides: {guides.shape}")

    required_cols = ["final_guide_id", "mlcb_gene", "mlcb_spacer", "final_gene_rank"]
    missing = [c for c in required_cols if c not in guides.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    genome_files = find_genome_fastas()

    if not genome_files:
        print("\nNo genome FASTA files found yet.")
        print(f"Put genome FASTA files here: {GENOME_DIR}")
        print("\nExpected extensions: .fa, .fasta, .fna")
        print("\nStill writing a guide-only preview file...")

    # Correct Cas-OFFinder query sequence:
    # 20-nt spacer + NNN, not spacer + NGG.
    guides["cas_offinder_query"] = guides["mlcb_spacer"].astype(str) + "NNN"
    guides["allowed_mismatches"] = ALLOWED_MISMATCHES

    preview_file = CONFIG_DIR / "cas_offinder_queries_preview.tsv"
    guides[
        [
            "final_guide_id",
            "mlcb_gene",
            "final_gene_rank",
            "mlcb_spacer",
            "cas_offinder_query",
            "allowed_mismatches",
            "functional_risk_penalized_score",
            "predicted_high_functional_offtarget_risk_probability",
        ]
    ].to_csv(preview_file, sep="\t", index=False)

    created_configs = []

    for genome_path in genome_files:
        genome_name = genome_path.stem
        config_file = CONFIG_DIR / f"cas_offinder_{genome_name}_top5_all4genes.txt"

        with open(config_file, "w") as f:
            f.write(str(genome_path) + "\n")
            f.write(CAS_OFFINDER_PATTERN + "\n")

            for _, row in guides.iterrows():
                query = row["cas_offinder_query"]
                f.write(f"{query} {ALLOWED_MISMATCHES}\n")

        created_configs.append(config_file)

    with open(SUMMARY_FILE, "w") as f:
        f.write("Cas-OFFinder Config File Summary\n")
        f.write("================================\n\n")
        f.write(f"Guide input file: {GUIDE_FILE}\n")
        f.write(f"Genome folder: {GENOME_DIR}\n")
        f.write(f"Cas-OFFinder pattern: {CAS_OFFINDER_PATTERN}\n")
        f.write(f"Allowed mismatches: {ALLOWED_MISMATCHES}\n\n")

        f.write(f"Total guides: {guides.shape[0]}\n")
        f.write(f"Genome FASTA files found: {len(genome_files)}\n\n")

        f.write("Guide query format used:\n")
        f.write("20-nt spacer + NNN\n\n")

        f.write("Created config files:\n")
        if created_configs:
            for cfg in created_configs:
                f.write(f"- {cfg}\n")
        else:
            f.write("- None yet, because no genome FASTA files were found.\n")

        f.write("\nSelected guides:\n")
        for _, row in guides.iterrows():
            f.write(
                f"{row['final_guide_id']} | {row['mlcb_gene']} | "
                f"rank={row['final_gene_rank']} | "
                f"query={row['cas_offinder_query']} | "
                f"mismatches={ALLOWED_MISMATCHES}\n"
            )

    print("\nWrote:")
    print(f"- {preview_file}")
    print(f"- {SUMMARY_FILE}")

    if created_configs:
        print("\nCreated config files:")
        for cfg in created_configs:
            print(f"- {cfg}")
    else:
        print("\nNo config files created yet because genome FASTA files are missing.")

    print("\nDone.")


if __name__ == "__main__":
    main()