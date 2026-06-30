#!/usr/bin/env python3

from pathlib import Path
import shutil

PROJECT_ROOT = Path(".").resolve()

DOWNLOAD_ROOT = PROJECT_ROOT / "data_cas_offinder/expanded_panel/ncbi_downloads/expanded_panel_36_accessions"
GENOME_OUT = PROJECT_ROOT / "data_cas_offinder/expanded_panel/genomes"
ANNOT_OUT = PROJECT_ROOT / "data_cas_offinder/expanded_panel/annotations"
META_OUT = PROJECT_ROOT / "data_cas_offinder/expanded_panel/metadata"

GENOME_OUT.mkdir(parents=True, exist_ok=True)
ANNOT_OUT.mkdir(parents=True, exist_ok=True)
META_OUT.mkdir(parents=True, exist_ok=True)

genome_files = sorted(DOWNLOAD_ROOT.rglob("*genomic.fna"))
gff_files = sorted(list(DOWNLOAD_ROOT.rglob("*.gff")) + list(DOWNLOAD_ROOT.rglob("*.gff3")))

def accession_from_path(path: Path):
    parts = path.parts
    for part in parts:
        if part.startswith("GCA_") or part.startswith("GCF_"):
            return part
    return path.stem.split("_genomic")[0]

copied_genomes = []
copied_gffs = []

for f in genome_files:
    acc = accession_from_path(f)
    out = GENOME_OUT / f"{acc}_genomic.fna"
    shutil.copy2(f, out)
    copied_genomes.append((acc, str(f), str(out)))

for f in gff_files:
    acc = accession_from_path(f)
    suffix = ".gff3" if f.name.endswith(".gff3") else ".gff"
    out = ANNOT_OUT / f"{acc}_genomic{suffix}"
    shutil.copy2(f, out)
    copied_gffs.append((acc, str(f), str(out)))

genome_accs = {x[0] for x in copied_genomes}
gff_accs = {x[0] for x in copied_gffs}

summary = []
summary.append("Expanded Panel Download Organization QC")
summary.append("=======================================")
summary.append(f"Download root: {DOWNLOAD_ROOT}")
summary.append(f"Genome output folder: {GENOME_OUT}")
summary.append(f"Annotation output folder: {ANNOT_OUT}")
summary.append("")
summary.append(f"Genome files copied: {len(copied_genomes)}")
summary.append(f"GFF/GFF3 files copied: {len(copied_gffs)}")
summary.append(f"Unique genome accessions: {len(genome_accs)}")
summary.append(f"Unique GFF accessions: {len(gff_accs)}")
summary.append(f"Accessions with genome but no GFF: {len(genome_accs - gff_accs)}")
summary.append(f"Accessions with GFF but no genome: {len(gff_accs - genome_accs)}")
summary.append("")
summary.append("QC status: " + ("PASS" if genome_accs == gff_accs and len(genome_accs) > 0 else "CHECK"))

summary_path = META_OUT / "expanded_panel_organization_qc.txt"
summary_path.write_text("\n".join(summary) + "\n")

print("\n".join(summary))
