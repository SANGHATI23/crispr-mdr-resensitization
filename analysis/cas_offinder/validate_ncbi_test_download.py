#!/usr/bin/env python3
"""
Validate one NCBI Datasets genome package before downloading the full panel.

Checks:
1. Exactly one genomic FASTA is present.
2. A GFF/GFF3 annotation file is present.
3. FASTA sequence IDs can be parsed.
4. GFF sequence IDs can be parsed.
5. Every GFF sequence ID occurs in the FASTA.
6. Basic genome and annotation counts are reported.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ASSEMBLY_ACCESSION = "GCA_000013465.1"

ASSEMBLY_DIR = (
    PROJECT_ROOT
    / "data_cas_offinder"
    / "expanded_panel"
    / "test_download"
    / "unzipped"
    / "ncbi_dataset"
    / "data"
    / ASSEMBLY_ACCESSION
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data_cas_offinder"
    / "expanded_panel"
    / "metadata"
    / "ncbi_test_download_validation.txt"
)


def fail(message: str) -> None:
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_fasta_ids(path: Path) -> list[str]:
    sequence_ids = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                sequence_id = line[1:].strip().split()[0]
                sequence_ids.append(sequence_id)

    return sequence_ids


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_id: str | None = None

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_id = line[1:].split()[0]
                lengths[current_id] = 0
            elif current_id is not None:
                lengths[current_id] += len(line)

    return lengths


def read_gff_ids_and_counts(
    path: Path,
) -> tuple[set[str], dict[str, int]]:
    sequence_ids: set[str] = set()
    feature_counts: dict[str, int] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            seqid = fields[0]
            feature_type = fields[2]

            sequence_ids.add(seqid)

            feature_counts[feature_type] = (
                feature_counts.get(feature_type, 0) + 1
            )

    return sequence_ids, feature_counts


def main() -> None:
    print("=" * 72)
    print("NCBI genome-package compatibility validation")
    print("=" * 72)
    print(f"Assembly directory: {ASSEMBLY_DIR}")

    if not ASSEMBLY_DIR.exists():
        fail(
            "Assembly directory not found. Confirm that the NCBI package "
            f"was unzipped correctly:\n{ASSEMBLY_DIR}"
        )

    fasta_files = sorted(
        ASSEMBLY_DIR.glob("*_genomic.fna")
    )

    gff_files = sorted(
        list(ASSEMBLY_DIR.glob("*.gff"))
        + list(ASSEMBLY_DIR.glob("*.gff3"))
    )

    if len(fasta_files) != 1:
        fail(
            "Expected exactly one genomic FASTA, but found "
            f"{len(fasta_files)}:\n{fasta_files}"
        )

    if len(gff_files) != 1:
        fail(
            "Expected exactly one GFF/GFF3 file, but found "
            f"{len(gff_files)}:\n{gff_files}"
        )

    fasta_file = fasta_files[0]
    gff_file = gff_files[0]

    fasta_ids = read_fasta_ids(fasta_file)
    fasta_lengths = read_fasta_lengths(fasta_file)

    gff_ids, feature_counts = read_gff_ids_and_counts(
        gff_file
    )

    if not fasta_ids:
        fail("No FASTA sequence identifiers were detected.")

    if not gff_ids:
        fail("No GFF sequence identifiers were detected.")

    fasta_id_set = set(fasta_ids)

    gff_ids_missing_from_fasta = sorted(
        gff_ids - fasta_id_set
    )

    fasta_ids_without_gff = sorted(
        fasta_id_set - gff_ids
    )

    compatible = len(gff_ids_missing_from_fasta) == 0

    total_sequence_length = sum(fasta_lengths.values())

    summary_lines = [
        "NCBI Test Download Validation",
        "=" * 72,
        "",
        f"Assembly accession: {ASSEMBLY_ACCESSION}",
        f"FASTA file: {fasta_file}",
        f"GFF file: {gff_file}",
        "",
        "FASTA statistics",
        "-" * 72,
        f"Number of FASTA records: {len(fasta_ids)}",
        f"Total sequence length: {total_sequence_length:,} bp",
        f"FASTA IDs: {sorted(fasta_id_set)}",
        "",
        "GFF statistics",
        "-" * 72,
        f"Number of GFF sequence IDs: {len(gff_ids)}",
        f"GFF sequence IDs: {sorted(gff_ids)}",
        f"Feature counts: {feature_counts}",
        "",
        "Compatibility checks",
        "-" * 72,
        (
            "GFF IDs absent from FASTA: "
            f"{gff_ids_missing_from_fasta}"
        ),
        (
            "FASTA IDs without GFF features: "
            f"{fasta_ids_without_gff}"
        ),
        (
            "FASTA/GFF identifier compatibility: "
            + ("PASS" if compatible else "FAIL")
        ),
        "",
        "Final status",
        "-" * 72,
        "PASS" if compatible else "FAIL",
    ]

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_FILE.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print(f"\nFASTA file: {fasta_file.name}")
    print(f"GFF file: {gff_file.name}")
    print(f"FASTA records: {len(fasta_ids)}")
    print(f"GFF sequence IDs: {len(gff_ids)}")
    print(f"Total sequence length: {total_sequence_length:,} bp")
    print(
        "FASTA/GFF compatibility: "
        + ("PASS" if compatible else "FAIL")
    )

    print(f"\nWrote: {OUTPUT_FILE.relative_to(PROJECT_ROOT)}")

    if not compatible:
        fail(
            "One or more GFF sequence identifiers were absent "
            "from the genomic FASTA."
        )


if __name__ == "__main__":
    main()