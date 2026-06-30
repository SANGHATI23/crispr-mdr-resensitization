#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HIT_LEVEL_CONTEXT = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/requirement3_hit_level_mobile_context_mapping.csv"

GENOME_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/genomes"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping/direct_plasmid_mapping"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_HIT = OUT_DIR / "requirement3_direct_plasmid_hit_level_mapping.csv"
OUT_GUIDE = OUT_DIR / "requirement3_direct_plasmid_guide_level_summary.csv"
OUT_SUMMARY = OUT_DIR / "requirement3_direct_plasmid_mapping_summary.txt"


PLASMID_TERMS = [
    "plasmid",
    "complete sequence plasmid",
    "circular plasmid",
]

CHROMOSOME_TERMS = [
    "chromosome",
    "complete genome",
    "chromosomal",
]


def extract_accession(path):
    m = re.search(r"(GC[AF]_\d+\.\d+)", path.name)
    if m:
        return m.group(1)
    return path.stem


def parse_fasta_headers(genome_dir):
    fasta_files = (
        list(genome_dir.glob("*.fna"))
        + list(genome_dir.glob("*.fa"))
        + list(genome_dir.glob("*.fasta"))
    )

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {genome_dir}")

    rows = []

    for fasta in fasta_files:
        accession = extract_accession(fasta)

        with open(fasta, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(">"):
                    continue

                header = line.strip().replace(">", "")
                seqid = header.split()[0]
                header_lower = header.lower()

                if any(term in header_lower for term in PLASMID_TERMS):
                    contig_class = "Plasmid_labeled_contig"
                    plasmid_flag = 1
                    chromosome_flag = 0
                elif any(term in header_lower for term in CHROMOSOME_TERMS):
                    contig_class = "Chromosome_labeled_contig"
                    plasmid_flag = 0
                    chromosome_flag = 1
                else:
                    contig_class = "Unclassified_contig"
                    plasmid_flag = 0
                    chromosome_flag = 0

                rows.append(
                    {
                        "assembly_accession": accession,
                        "seqid": seqid,
                        "fasta_header": header,
                        "direct_contig_class": contig_class,
                        "direct_plasmid_contig_flag": plasmid_flag,
                        "direct_chromosome_contig_flag": chromosome_flag,
                        "genome_fasta_file": str(fasta.relative_to(PROJECT_ROOT)),
                    }
                )

    return pd.DataFrame(rows)


def pick_col(df, candidates, required=True):
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    if required:
        raise ValueError(
            f"Missing required column. Tried {candidates}. Available: {list(df.columns)}"
        )
    return None


def main():
    print("Requirement 3 direct plasmid/chromosome contig mapping started")
    print(f"Project root: {PROJECT_ROOT}")

    if not HIT_LEVEL_CONTEXT.exists():
        raise FileNotFoundError(f"Hit-level context table not found: {HIT_LEVEL_CONTEXT}")

    hits = pd.read_csv(HIT_LEVEL_CONTEXT)
    print(f"Loaded hit-level Requirement 3 table: {hits.shape}")

    headers = parse_fasta_headers(GENOME_DIR)
    print(f"Parsed FASTA headers: {headers.shape}")

    accession_col = pick_col(
        hits,
        [
            "requirement3_accession_used",
            "accession",
            "genome_accession",
            "assembly_accession",
        ],
        required=True,
    )

    seqid_col = pick_col(
        hits,
        [
            "requirement3_seqid_used",
            "seqid_for_mapping",
            "seqid",
            "contig",
            "contig_id",
        ],
        required=True,
    )

    guide_col = pick_col(
        hits,
        ["guide_id", "selected_guide_id", "guide"],
        required=True,
    )

    gene_col = pick_col(
        hits,
        ["target_gene", "gene", "amr_gene"],
        required=True,
    )

    hits = hits.rename(
        columns={
            accession_col: "assembly_accession",
            seqid_col: "seqid",
            guide_col: "guide_id",
            gene_col: "target_gene",
        }
    )

    merged = hits.merge(
        headers,
        on=["assembly_accession", "seqid"],
        how="left",
    )

    merged["direct_contig_class"] = merged["direct_contig_class"].fillna("No_FASTA_header_match")
    merged["direct_plasmid_contig_flag"] = merged["direct_plasmid_contig_flag"].fillna(0).astype(int)
    merged["direct_chromosome_contig_flag"] = merged["direct_chromosome_contig_flag"].fillna(0).astype(int)

    merged.to_csv(OUT_HIT, index=False)

    guide_summary = (
        merged.groupby(["target_gene", "guide_id"], dropna=False)
        .agg(
            total_hit_rows=("guide_id", "count"),
            direct_plasmid_contig_hits=("direct_plasmid_contig_flag", "sum"),
            direct_chromosome_contig_hits=("direct_chromosome_contig_flag", "sum"),
            no_fasta_header_match_hits=("direct_contig_class", lambda x: int((x == "No_FASTA_header_match").sum())),
            unclassified_contig_hits=("direct_contig_class", lambda x: int((x == "Unclassified_contig").sum())),
        )
        .reset_index()
    )

    guide_summary["direct_plasmid_hit_fraction"] = (
        guide_summary["direct_plasmid_contig_hits"] / guide_summary["total_hit_rows"]
    ).fillna(0)

    def classify(row):
        if row["direct_plasmid_contig_hits"] > 0:
            return "Direct_plasmid_contig_evidence_present"
        if row["no_fasta_header_match_hits"] == row["total_hit_rows"]:
            return "No_direct_contig_header_match"
        return "No_direct_plasmid_contig_evidence_detected"

    guide_summary["direct_plasmid_mapping_decision"] = guide_summary.apply(classify, axis=1)

    guide_summary.to_csv(OUT_GUIDE, index=False)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Requirement 3 Direct Plasmid/Chromosome Contig Mapping Summary\n")
        f.write("=============================================================\n\n")

        f.write(f"Input hit-level context table: {HIT_LEVEL_CONTEXT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Genome FASTA folder: {GENOME_DIR.relative_to(PROJECT_ROOT)}\n\n")

        f.write(f"Hit rows evaluated: {len(merged)}\n")
        f.write(f"FASTA headers parsed: {len(headers)}\n")
        f.write(f"Guide rows summarized: {len(guide_summary)}\n\n")

        f.write("Hit-level direct contig class counts:\n")
        for k, v in merged["direct_contig_class"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nGuide-level direct plasmid mapping decision counts:\n")
        for k, v in guide_summary["direct_plasmid_mapping_decision"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nGuides with direct plasmid-contig evidence:\n")
        subset = guide_summary[guide_summary["direct_plasmid_contig_hits"] > 0]
        if subset.empty:
            f.write("- None detected from FASTA header labels.\n")
        else:
            for _, r in subset.iterrows():
                f.write(
                    f"- {r['target_gene']} | {r['guide_id']} | "
                    f"plasmid_contig_hits={int(r['direct_plasmid_contig_hits'])} | "
                    f"total_hits={int(r['total_hit_rows'])} | "
                    f"fraction={r['direct_plasmid_hit_fraction']:.3f}\n"
                )

        f.write("\nOutput files:\n")
        f.write(f"- {OUT_HIT.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_GUIDE.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {OUT_SUMMARY.relative_to(PROJECT_ROOT)}\n")

    print(f"Wrote: {OUT_HIT.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_GUIDE.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {OUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nDirect contig class counts:")
    print(merged["direct_contig_class"].value_counts().to_string())

    print("\nGuides with direct plasmid-contig evidence:")
    subset = guide_summary[guide_summary["direct_plasmid_contig_hits"] > 0]
    if subset.empty:
        print("None detected from FASTA header labels.")
    else:
        print(subset.to_string(index=False))

    print("\nRequirement 3 direct plasmid/chromosome contig mapping completed.")


if __name__ == "__main__":
    main()
