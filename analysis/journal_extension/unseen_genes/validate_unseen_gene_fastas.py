from pathlib import Path
import pandas as pd
import re

ROOT = Path(__file__).resolve().parents[3]

REFERENCE_PLAN = ROOT / "results_journal_extension" / "unseen_genes" / "unseen_gene_reference_plan.csv"
FASTA_DIR = ROOT / "data_journal_extension" / "unseen_genes"

OUT_DIR = ROOT / "results_journal_extension" / "unseen_genes"
OUT_TABLE = OUT_DIR / "unseen_gene_fasta_validation.csv"
OUT_SUMMARY = OUT_DIR / "unseen_gene_fasta_validation_summary.txt"

VALID_BASES = set("ACGTN")


def read_fasta(path):
    if not path.exists():
        return None, ""

    header = None
    seq_parts = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
                else:
                    # Only first record is used for now.
                    pass
            else:
                seq_parts.append(line.upper())

    seq = "".join(seq_parts)
    seq = re.sub(r"\s+", "", seq)

    return header, seq


def validate_sequence(seq):
    if not seq:
        return {
            "sequence_present": "No",
            "length_bp": 0,
            "valid_dna_alphabet": "No",
            "length_divisible_by_3": "No",
            "starts_with_atg": "No",
            "has_terminal_stop": "No",
            "validation_status": "Missing_sequence",
            "validation_notes": "No sequence found."
        }

    invalid_chars = sorted(set(seq) - VALID_BASES)
    valid_alpha = len(invalid_chars) == 0
    length_bp = len(seq)
    div3 = length_bp % 3 == 0
    starts_atg = seq.startswith("ATG")
    terminal_stop = seq[-3:] in {"TAA", "TAG", "TGA"} if length_bp >= 3 else False

    notes = []
    status = "Pass_with_warnings"

    if not valid_alpha:
        status = "Fail"
        notes.append(f"Invalid characters found: {''.join(invalid_chars)}")

    if length_bp < 100:
        status = "Fail"
        notes.append("Sequence length is unexpectedly short for AMR coding sequence.")

    if not div3:
        notes.append("Length is not divisible by 3; verify CDS boundaries.")

    if not starts_atg:
        notes.append("Sequence does not start with ATG; verify CDS start.")

    if not terminal_stop:
        notes.append("Sequence does not end with a standard stop codon; verify CDS end.")

    if valid_alpha and length_bp >= 100 and div3 and starts_atg and terminal_stop:
        status = "Pass"

    return {
        "sequence_present": "Yes",
        "length_bp": length_bp,
        "valid_dna_alphabet": "Yes" if valid_alpha else "No",
        "length_divisible_by_3": "Yes" if div3 else "No",
        "starts_with_atg": "Yes" if starts_atg else "No",
        "has_terminal_stop": "Yes" if terminal_stop else "No",
        "validation_status": status,
        "validation_notes": "; ".join(notes) if notes else "No issues detected."
    }


def main():
    plan = pd.read_csv(REFERENCE_PLAN)

    rows = []

    for _, r in plan.iterrows():
        gene = r["gene"]
        fasta_path = ROOT / r["planned_fasta_path"]

        header, seq = read_fasta(fasta_path)
        val = validate_sequence(seq)

        rows.append({
            "gene": gene,
            "fasta_path": str(fasta_path.relative_to(ROOT)),
            "file_exists": "Yes" if fasta_path.exists() else "No",
            "header": header if header else "",
            **val
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TABLE, index=False)

    lines = []
    lines.append("Unseen Gene FASTA Validation Summary")
    lines.append("====================================")
    lines.append("")
    lines.append(f"Reference plan: {REFERENCE_PLAN}")
    lines.append(f"FASTA directory: {FASTA_DIR}")
    lines.append(f"Output table: {OUT_TABLE}")
    lines.append("")
    lines.append("Validation status counts:")
    for k, v in df["validation_status"].value_counts(dropna=False).items():
        lines.append(f"- {k}: {v}")

    lines.append("")
    lines.append("Per-gene status:")
    for _, r in df.iterrows():
        lines.append(
            f"- {r['gene']}: file_exists={r['file_exists']} | "
            f"length={r['length_bp']} | status={r['validation_status']} | "
            f"notes={r['validation_notes']}"
        )

    OUT_SUMMARY.write_text("\n".join(lines))

    print(f"Wrote: {OUT_TABLE}")
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
