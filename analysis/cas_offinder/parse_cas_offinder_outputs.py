import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "results_cas_offinder" / "raw"
PARSED_DIR = PROJECT_ROOT / "results_cas_offinder" / "parsed"
PARSED_DIR.mkdir(parents=True, exist_ok=True)

GUIDE_FILE = PROJECT_ROOT / "results_final_selection" / "top5_per_gene_risk_aware_guides.csv"

OUTPUT_FILE = PARSED_DIR / "cas_offinder_all_genome_offtarget_hits_parsed.csv"
SUMMARY_FILE = PARSED_DIR / "cas_offinder_parsed_summary.txt"

RAW_FILES = {
    "Klebsiella_pneumoniae_reference": RAW_DIR / "klebsiella_top5_all4genes_offtargets.txt",
    "Ecoli_reference": RAW_DIR / "ecoli_top5_all4genes_offtargets.txt",
    "Staphylococcus_aureus_MRSA_reference": RAW_DIR / "staph_aureus_top5_all4genes_offtargets.txt",
}


def parse_line(line, genome_label):
    """
    Cas-OFFinder output format approximately:
    query  chromosome/contig description...  position  off_target_sequence  strand  mismatch_count

    Because the contig description contains spaces, we parse from both ends:
    first token = query
    last token = mismatch_count
    second-last token = strand
    third-last token = off_target_sequence
    fourth-last token = position
    everything between query and position = chromosome/contig description
    """
    parts = line.strip().split()

    if len(parts) < 6:
        return None

    query = parts[0]
    mismatch_count = parts[-1]
    strand = parts[-2]
    off_target_sequence = parts[-3]
    position = parts[-4]
    contig_description = " ".join(parts[1:-4])

    try:
        mismatch_count = int(mismatch_count)
    except ValueError:
        mismatch_count = None

    try:
        position = int(position)
    except ValueError:
        position = None

    spacer = query[:20]

    return {
        "genome_label": genome_label,
        "cas_offinder_query": query,
        "query_spacer": spacer,
        "contig_description": contig_description,
        "position": position,
        "off_target_sequence": off_target_sequence,
        "strand": strand,
        "mismatch_count": mismatch_count,
    }


def main():
    print("Parsing Cas-OFFinder raw output files...")

    if not GUIDE_FILE.exists():
        raise FileNotFoundError(f"Missing guide file: {GUIDE_FILE}")

    guides = pd.read_csv(GUIDE_FILE)
    print(f"Loaded guide table: {guides.shape}")

    required_guide_cols = [
        "final_guide_id",
        "mlcb_gene",
        "mlcb_spacer",
        "final_gene_rank",
        "functional_risk_penalized_score",
        "predicted_high_functional_offtarget_risk_probability",
    ]
    missing = [c for c in required_guide_cols if c not in guides.columns]
    if missing:
        raise ValueError(f"Guide table missing columns: {missing}")

    guide_map = guides[required_guide_cols].copy()
    guide_map = guide_map.rename(columns={"mlcb_spacer": "query_spacer"})

    rows = []

    for genome_label, raw_file in RAW_FILES.items():
        if not raw_file.exists():
            print(f"WARNING: missing raw file: {raw_file}")
            continue

        print(f"Reading: {raw_file}")

        with open(raw_file, "r") as f:
            for line in f:
                if not line.strip():
                    continue

                parsed = parse_line(line, genome_label)
                if parsed is not None:
                    rows.append(parsed)

    if not rows:
        raise ValueError("No Cas-OFFinder hits parsed. Check raw output files.")

    hits = pd.DataFrame(rows)
    print(f"Parsed raw hits: {hits.shape}")

    # Merge hits back to selected guides
    merged = hits.merge(guide_map, on="query_spacer", how="left")

    # Reorder columns
    cols = [
        "final_guide_id",
        "mlcb_gene",
        "final_gene_rank",
        "query_spacer",
        "cas_offinder_query",
        "genome_label",
        "contig_description",
        "position",
        "strand",
        "off_target_sequence",
        "mismatch_count",
        "functional_risk_penalized_score",
        "predicted_high_functional_offtarget_risk_probability",
    ]

    merged = merged[cols].copy()

    # Flag unmatched rows
    merged["guide_mapping_status"] = merged["final_guide_id"].apply(
        lambda x: "Mapped_to_selected_guide" if pd.notna(x) else "Unmapped_check_spacer"
    )

    # Sort for readability
    merged = merged.sort_values(
        by=["mlcb_gene", "final_gene_rank", "genome_label", "mismatch_count", "position"],
        na_position="last"
    )

    merged.to_csv(OUTPUT_FILE, index=False)

    # Summary tables
    total_hits = len(merged)
    hits_by_genome = merged.groupby("genome_label").size().reset_index(name="hit_count")
    hits_by_gene = merged.groupby("mlcb_gene").size().reset_index(name="hit_count")
    hits_by_mismatch = merged.groupby("mismatch_count").size().reset_index(name="hit_count")
    hits_by_guide = (
        merged.groupby(["final_guide_id", "mlcb_gene", "final_gene_rank"])
        .size()
        .reset_index(name="hit_count")
        .sort_values(["mlcb_gene", "final_gene_rank"])
    )

    with open(SUMMARY_FILE, "w") as f:
        f.write("Parsed Cas-OFFinder Output Summary\n")
        f.write("=================================\n\n")
        f.write(f"Total parsed off-target hits: {total_hits}\n")
        f.write(f"Output CSV: {OUTPUT_FILE}\n\n")

        f.write("Hits by genome:\n")
        f.write(hits_by_genome.to_string(index=False))
        f.write("\n\n")

        f.write("Hits by gene:\n")
        f.write(hits_by_gene.to_string(index=False))
        f.write("\n\n")

        f.write("Hits by mismatch count:\n")
        f.write(hits_by_mismatch.to_string(index=False))
        f.write("\n\n")

        f.write("Hits by guide:\n")
        f.write(hits_by_guide.to_string(index=False))
        f.write("\n\n")

        unmapped = merged[merged["guide_mapping_status"] != "Mapped_to_selected_guide"]
        f.write(f"Unmapped rows: {len(unmapped)}\n")

    print("\nWrote:")
    print(f"- {OUTPUT_FILE}")
    print(f"- {SUMMARY_FILE}")

    print("\nHits by genome:")
    print(hits_by_genome.to_string(index=False))

    print("\nHits by gene:")
    print(hits_by_gene.to_string(index=False))

    print("\nHits by mismatch count:")
    print(hits_by_mismatch.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()