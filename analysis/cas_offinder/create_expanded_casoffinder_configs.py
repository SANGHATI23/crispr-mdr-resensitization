#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(".").resolve()

GUIDE_TSV = PROJECT_ROOT / "results_final_selection/cas_offinder_input_guides_top5_per_gene.tsv"
GUIDE_CSV = PROJECT_ROOT / "results_final_selection/top5_per_gene_risk_aware_guides.csv"

GENOME_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/genomes"
CONFIG_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/configs"
RAW_OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/raw"
META_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/metadata"

CONFIG_DIR.mkdir(parents=True, exist_ok=True)
RAW_OUT_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

SPCAS9_PATTERN = "NNNNNNNNNNNNNNNNNNNNNGG"
MISMATCHES = 4


def load_guides():
    if GUIDE_TSV.exists():
        df = pd.read_csv(GUIDE_TSV, sep="\t")
    elif GUIDE_CSV.exists():
        df = pd.read_csv(GUIDE_CSV)
    else:
        raise FileNotFoundError("Could not find guide TSV or guide CSV in results_final_selection/")

    cols = {c.lower(): c for c in df.columns}

    spacer_col = None
    for possible in ["spacer", "mlcb_spacer", "guide", "guide_sequence", "sequence"]:
        if possible in cols:
            spacer_col = cols[possible]
            break

    if spacer_col is None:
        raise ValueError(f"Could not detect spacer column. Columns found: {list(df.columns)}")

    guide_id_col = None
    for possible in ["guide_id", "final_guide_id", "id", "candidate_id"]:
        if possible in cols:
            guide_id_col = cols[possible]
            break

    gene_col = None
    for possible in ["gene", "mlcb_gene", "target_gene", "amr_gene"]:
        if possible in cols:
            gene_col = cols[possible]
            break

    guides = df.copy()
    guides["spacer_clean"] = guides[spacer_col].astype(str).str.upper().str.strip()

    guides = guides[guides["spacer_clean"].str.len() == 20].copy()
    guides = guides.drop_duplicates(subset=["spacer_clean"]).copy()

    if guide_id_col is None:
        guides["guide_id_clean"] = ["guide_" + str(i + 1) for i in range(len(guides))]
    else:
        guides["guide_id_clean"] = guides[guide_id_col].astype(str)

    if gene_col is None:
        guides["gene_clean"] = "unknown"
    else:
        guides["gene_clean"] = guides[gene_col].astype(str)

    guides["casoffinder_query"] = guides["spacer_clean"] + "NNN"

    bad = guides[guides["casoffinder_query"].str.len() != 23]
    if len(bad) > 0:
        raise ValueError("Some Cas-OFFinder queries are not 23 nt after spacer+NNN formatting.")

    return guides[["gene_clean", "guide_id_clean", "spacer_clean", "casoffinder_query"]]


def accession_from_genome_file(path: Path):
    name = path.name
    if name.endswith("_genomic.fna"):
        return name.replace("_genomic.fna", "")
    return path.stem


def main():
    guides = load_guides()
    genome_files = sorted(GENOME_DIR.glob("*_genomic.fna"))

    if len(guides) == 0:
        raise ValueError("No guides loaded.")
    if len(genome_files) == 0:
        raise ValueError("No genome FASTA files found.")

    manifest_rows = []

    for genome_file in genome_files:
        acc = accession_from_genome_file(genome_file)
        config_file = CONFIG_DIR / f"{acc}_casoffinder_config.txt"
        raw_output_file = RAW_OUT_DIR / f"{acc}_expanded_offtargets.txt"

        lines = []
        lines.append(str(genome_file))
        lines.append(SPCAS9_PATTERN)

        for _, row in guides.iterrows():
            lines.append(f"{row['casoffinder_query']} {MISMATCHES}")

        config_file.write_text("\n".join(lines) + "\n")

        manifest_rows.append({
            "accession": acc,
            "genome_file": str(genome_file),
            "config_file": str(config_file),
            "raw_output_file": str(raw_output_file),
            "guide_count": len(guides),
            "pattern": SPCAS9_PATTERN,
            "mismatches": MISMATCHES,
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = META_DIR / "expanded_casoffinder_run_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    run_script = PROJECT_ROOT / "analysis/cas_offinder/run_expanded_casoffinder_all.sh"
    script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "mkdir -p results_cas_offinder/expanded_panel/raw",
        "",
    ]

    for _, row in manifest.iterrows():
        script_lines.append(f"echo 'Running Cas-OFFinder for {row['accession']}'")
        script_lines.append(f"cas-offinder {row['config_file']} C {row['raw_output_file']}")
        script_lines.append("")

    run_script.write_text("\n".join(script_lines) + "\n")

    qc_path = META_DIR / "expanded_casoffinder_config_qc.txt"

    summary = []
    summary.append("Expanded Cas-OFFinder Config QC")
    summary.append("===============================")
    summary.append(f"Guides loaded: {len(guides)}")
    summary.append(f"Genome FASTA files found: {len(genome_files)}")
    summary.append(f"Config files created: {len(manifest)}")
    summary.append(f"Mismatch threshold: {MISMATCHES}")
    summary.append(f"SpCas9 pattern: {SPCAS9_PATTERN}")
    summary.append(f"Manifest: {manifest_path}")
    summary.append(f"Run script: {run_script}")
    summary.append("")
    summary.append("Guide counts by gene:")
    summary.append(str(guides["gene_clean"].value_counts(dropna=False)))
    summary.append("")
    summary.append("QC status: " + ("PASS" if len(guides) == 20 and len(manifest) == len(genome_files) else "CHECK"))

    qc_path.write_text("\n".join(summary) + "\n")

    print("\n".join(summary))


if __name__ == "__main__":
    main()
