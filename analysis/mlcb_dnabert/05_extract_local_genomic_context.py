#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def clean(seq):
    if pd.isna(seq):
        return ""
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("\r", "")


def rc(seq):
    return str(Seq(seq).reverse_complement())


def pick_col(df, options):
    for c in options:
        if c in df.columns:
            return c
    return None


def read_fastas(fasta_files):
    records = []
    for fp in fasta_files:
        p = Path(fp)
        if not p.exists():
            print(f"WARNING missing FASTA: {fp}")
            continue
        for rec in SeqIO.parse(str(p), "fasta"):
            records.append({
                "source_file": str(p),
                "record_id": rec.id,
                "description": rec.description,
                "sequence": clean(str(rec.seq))
            })
    return records


def extract_window(seq, start, end, context_bp):
    mid = (start + end) // 2
    left = max(0, mid - context_bp // 2)
    right = min(len(seq), left + context_bp)

    if right - left < context_bp:
        left = max(0, right - context_bp)

    return seq[left:right], left, right


def gene_aliases(gene):
    g = str(gene).lower().replace("-", "").replace("_", "").replace(" ", "")

    alias_map = {
        "blakpc": ["blakpc", "kpc"],
        "blandm1": ["blandm", "blandm1", "ndm", "ndm1"],
        "blandm": ["blandm", "ndm"],
        "mcr1": ["mcr1", "mcr-1"],
        "meca": ["meca"],
    }

    return alias_map.get(g, [g])


def record_matches_gene(record, gene):
    text = " ".join([
        record["source_file"],
        record["record_id"],
        record["description"]
    ]).lower()

    return any(a in text for a in gene_aliases(gene))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results_mlcb/tables/mlcb_guide_dataset.csv")
    parser.add_argument("--output", default="results_mlcb/tables/mlcb_guide_dataset_with_local_context.csv")
    parser.add_argument("--context-bp", type=int, default=100)
    args = parser.parse_args()

    fasta_files = [
        "data/targets_multistrain/blaKPC_multi.fasta",
        "data/targets_multistrain/blaNDM1_multi.fasta",
        "data/targets_multistrain/mcr1_multi.fasta",
        "data/targets_multistrain/mecA_multi.fasta",
        "data/plasmids/blaKPC_plasmid.fasta",
        "data/plasmids/blaNDM_plasmid.fasta",
        "data/plasmids/mcr1_plasmid.fasta",
        "data/genomes/staphylococcus_aureus_genome.fasta",
        "data/genomes/ecoli_genome.fasta",
    ]

    df = pd.read_csv(args.input)

    spacer_col = pick_col(df, ["mlcb_spacer", "spacer", "guide_sequence", "guide", "sgRNA", "sequence"])
    gene_col = pick_col(df, ["mlcb_gene", "gene", "target_gene", "Gene"])

    if spacer_col is None:
        raise ValueError("No spacer column found.")

    records = read_fastas(fasta_files)

    print("Loaded dataset:", df.shape)
    print("Spacer column:", spacer_col)
    print("Gene column:", gene_col)
    print("Reference records:", len(records))

    if len(records) == 0:
        raise ValueError("No FASTA records loaded.")

    output_rows = []

    for _, row in df.iterrows():
        spacer = clean(row[spacer_col])
        gene = row[gene_col] if gene_col else ""

        best = None
        all_matches = []

        if len(spacer) >= 15:
            spacer_rc = rc(spacer)

            for rec in records:
                seq = rec["sequence"]

                pos = seq.find(spacer)
                if pos != -1:
                    all_matches.append((rec, pos, pos + len(spacer), "+"))

                pos_rc = seq.find(spacer_rc)
                if pos_rc != -1:
                    all_matches.append((rec, pos_rc, pos_rc + len(spacer), "-"))

            gene_matches = [m for m in all_matches if record_matches_gene(m[0], gene)]
            best = gene_matches[0] if gene_matches else (all_matches[0] if all_matches else None)

        if best is None:
            output_rows.append({
                "local_context_sequence": np.nan,
                "local_context_length": 0,
                "context_status": "not_found_in_reference",
                "target_record_id": np.nan,
                "target_source_file": np.nan,
                "target_start_0based": np.nan,
                "target_end_0based": np.nan,
                "target_strand": np.nan,
                "context_start_0based": np.nan,
                "context_end_0based": np.nan,
                "number_of_reference_matches": len(all_matches),
            })
        else:
            rec, start, end, strand = best
            context, cstart, cend = extract_window(rec["sequence"], start, end, args.context_bp)

            if strand == "-":
                context = rc(context)

            output_rows.append({
                "local_context_sequence": context,
                "local_context_length": len(context),
                "context_status": "found",
                "target_record_id": rec["record_id"],
                "target_source_file": rec["source_file"],
                "target_start_0based": start,
                "target_end_0based": end,
                "target_strand": strand,
                "context_start_0based": cstart,
                "context_end_0based": cend,
                "number_of_reference_matches": len(all_matches),
            })

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(output_rows)], axis=1)

    if "dnabert_sequence" in out.columns:
        out["dnabert_sequence_original"] = out["dnabert_sequence"]
    else:
        out["dnabert_sequence_original"] = out[spacer_col]

    out["dnabert_sequence"] = out["local_context_sequence"].fillna(out[spacer_col].astype(str))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    summary_path = "results_mlcb/tables/local_context_extraction_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Local genomic context extraction summary\n")
        f.write(f"Input rows: {len(out)}\n")
        f.write(str(out["context_status"].value_counts()))

    print("Done.")
    print("Wrote:", args.output)
    print("Wrote:", summary_path)
    print(out["context_status"].value_counts())


if __name__ == "__main__":
    main()
