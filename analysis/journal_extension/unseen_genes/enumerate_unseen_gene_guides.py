from pathlib import Path
import pandas as pd
import re

ROOT = Path(__file__).resolve().parents[3]

FASTA_DIR = ROOT / "data_journal_extension" / "unseen_genes"
OUT_DIR = ROOT / "results_journal_extension" / "unseen_genes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "unseen_gene_candidate_guides_all.csv"
OUT_SUMMARY = OUT_DIR / "unseen_gene_candidate_guides_summary.txt"

GENES = ["vanA", "tetM", "ermB", "aac6Ib", "qnrS"]

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    header = ""
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
            else:
                seq_parts.append(line.upper())
    return header, "".join(seq_parts)


def revcomp(seq):
    return seq.translate(COMP)[::-1].upper()


def gc_content(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def simple_on_target_score(spacer):
    """
    Conservative placeholder-style on-target proxy for unseen-gene extension.
    This is not a replacement for RS3/Azimuth; it gives an initial candidate screen.
    """
    gc = gc_content(spacer)

    # Best around 40-60 GC; penalize extremes.
    gc_score = max(0, 100 - abs(gc - 50) * 2)

    # Avoid poly-T and long homopolymers.
    homopolymer_penalty = 0
    for base in "ACGT":
        if base * 4 in spacer:
            homopolymer_penalty += 10
        if base * 5 in spacer:
            homopolymer_penalty += 20

    score = max(0, min(100, gc_score - homopolymer_penalty))
    return round(score, 2)


def classify_gc(gc):
    if gc < 30:
        return "Low_GC"
    if gc > 70:
        return "High_GC"
    return "Acceptable_GC"


def target_region(position, length):
    frac = position / max(length, 1)
    if frac <= 0.33:
        return "early_coding_region"
    if frac <= 0.66:
        return "middle_coding_region"
    return "late_coding_region"


def enumerate_plus_strand(seq, gene, header):
    rows = []
    # Pattern: N20 + NGG, PAM positions i+20:i+23
    for i in range(0, len(seq) - 23 + 1):
        spacer = seq[i:i+20]
        pam = seq[i+20:i+23]
        if len(spacer) == 20 and len(pam) == 3 and pam[1:] == "GG" and "N" not in spacer + pam:
            gc = gc_content(spacer)
            rows.append({
                "gene": gene,
                "source_header": header,
                "strand": "+",
                "position_1based": i + 1,
                "spacer": spacer,
                "pam": pam,
                "protospacer_plus_pam": spacer + pam,
                "gc_content": round(gc, 2),
                "gc_class": classify_gc(gc),
                "target_region_bin": target_region(i + 1, len(seq)),
                "simple_on_target_proxy": simple_on_target_score(spacer),
                "sequence_length_bp": len(seq),
            })
    return rows


def enumerate_minus_strand(seq, gene, header):
    rows = []
    rc = revcomp(seq)

    # Enumerate guides on reverse complement, then map approximate position back.
    for i in range(0, len(rc) - 23 + 1):
        spacer = rc[i:i+20]
        pam = rc[i+20:i+23]
        if len(spacer) == 20 and len(pam) == 3 and pam[1:] == "GG" and "N" not in spacer + pam:
            # Approximate genomic coordinate on original plus strand
            original_start_1based = len(seq) - (i + 23) + 1
            gc = gc_content(spacer)
            rows.append({
                "gene": gene,
                "source_header": header,
                "strand": "-",
                "position_1based": original_start_1based,
                "spacer": spacer,
                "pam": pam,
                "protospacer_plus_pam": spacer + pam,
                "gc_content": round(gc, 2),
                "gc_class": classify_gc(gc),
                "target_region_bin": target_region(original_start_1based, len(seq)),
                "simple_on_target_proxy": simple_on_target_score(spacer),
                "sequence_length_bp": len(seq),
            })
    return rows


def main():
    all_rows = []

    for gene in GENES:
        fasta = FASTA_DIR / f"{gene}.fa"
        if not fasta.exists():
            print(f"Skipping {gene}: missing {fasta}")
            continue

        header, seq = read_fasta(fasta)
        rows = enumerate_plus_strand(seq, gene, header)
        rows += enumerate_minus_strand(seq, gene, header)

        gene_df = pd.DataFrame(rows)
        gene_out = OUT_DIR / f"{gene}_candidate_guides.csv"
        gene_df.to_csv(gene_out, index=False)
        print(f"Wrote: {gene_out} ({len(gene_df)} guides)")

        all_rows.extend(rows)

    all_df = pd.DataFrame(all_rows)
    all_df.to_csv(OUT_ALL, index=False)

    lines = []
    lines.append("Unseen Gene Candidate Guide Enumeration Summary")
    lines.append("==============================================")
    lines.append("")
    lines.append(f"Output all-guides table: {OUT_ALL}")
    lines.append(f"Total candidate guides: {len(all_df)}")
    lines.append("")

    if len(all_df) > 0:
        lines.append("Guides by gene:")
        for gene, count in all_df["gene"].value_counts().sort_index().items():
            lines.append(f"- {gene}: {count}")

        lines.append("")
        lines.append("Guides by gene and strand:")
        strand_counts = all_df.groupby(["gene", "strand"]).size().reset_index(name="count")
        for _, r in strand_counts.iterrows():
            lines.append(f"- {r['gene']} strand {r['strand']}: {r['count']}")

        lines.append("")
        lines.append("Top 5 guides by simple on-target proxy:")
        top = all_df.sort_values("simple_on_target_proxy", ascending=False).head(5)
        for _, r in top.iterrows():
            lines.append(
                f"- {r['gene']} | {r['spacer']} | PAM={r['pam']} | "
                f"strand={r['strand']} | pos={r['position_1based']} | "
                f"score={r['simple_on_target_proxy']} | region={r['target_region_bin']}"
            )

    lines.append("")
    lines.append("Interpretation:")
    lines.append(
        "This script creates an initial SpCas9 NGG candidate-guide set for prospective unseen AMR genes. "
        "The simple_on_target_proxy is only a first-pass sequence suitability proxy and should be replaced "
        "or supplemented by RS3/Azimuth-style scoring before manuscript-level claims."
    )

    OUT_SUMMARY.write_text("\n".join(lines))
    print(f"Wrote: {OUT_SUMMARY}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
