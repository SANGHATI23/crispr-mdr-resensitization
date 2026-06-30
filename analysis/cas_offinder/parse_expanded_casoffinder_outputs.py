#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(".").resolve()

RAW_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/raw"
GUIDE_FILE = PROJECT_ROOT / "results_final_selection/cas_offinder_input_guides_top5_per_gene.tsv"
MANIFEST_FILE = PROJECT_ROOT / "data_cas_offinder/expanded_panel/metadata/expanded_casoffinder_run_manifest.csv"

OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/parsed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_HITS = OUT_DIR / "expanded_casoffinder_all_genome_offtarget_hits_parsed.csv"
OUT_SUMMARY = OUT_DIR / "expanded_casoffinder_parsed_summary.txt"

def load_guides():
    guides = pd.read_csv(GUIDE_FILE, sep="\t")
    guides.columns = [c.strip() for c in guides.columns]

    if "mlcb_spacer" in guides.columns:
        spacer_col = "mlcb_spacer"
    elif "spacer" in guides.columns:
        spacer_col = "spacer"
    else:
        raise ValueError(f"Cannot find spacer column. Columns: {list(guides.columns)}")

    if "final_guide_id" in guides.columns:
        guide_id_col = "final_guide_id"
    elif "guide_id" in guides.columns:
        guide_id_col = "guide_id"
    else:
        guide_id_col = None

    if "mlcb_gene" in guides.columns:
        gene_col = "mlcb_gene"
    elif "gene" in guides.columns:
        gene_col = "gene"
    else:
        gene_col = None

    guides["spacer"] = guides[spacer_col].astype(str).str.upper().str.strip()
    guides["casoffinder_query"] = guides["spacer"] + "NNN"

    if guide_id_col:
        guides["guide_id"] = guides[guide_id_col].astype(str)
    else:
        guides["guide_id"] = ["guide_" + str(i + 1) for i in range(len(guides))]

    if gene_col:
        guides["target_gene"] = guides[gene_col].astype(str)
    else:
        guides["target_gene"] = "unknown"

    return guides[["guide_id", "target_gene", "spacer", "casoffinder_query"]].drop_duplicates()

def parse_raw_file(path: Path, accession: str):
    rows = []

    with path.open() as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 6:
                rows.append({
                    "accession": accession,
                    "source_file": str(path),
                    "line_num": line_num,
                    "parse_status": "FAILED_TOO_FEW_FIELDS",
                    "raw_line": line,
                })
                continue

            try:
                query_sequence = parts[0]
                mismatch_count = int(parts[-1])
                strand = parts[-2]
                off_target_sequence = parts[-3]
                position_1based = int(parts[-4])
                contig_description = " ".join(parts[1:-4])

                rows.append({
                    "accession": accession,
                    "source_file": str(path),
                    "line_num": line_num,
                    "parse_status": "PASS",
                    "query_sequence": query_sequence,
                    "contig_description": contig_description,
                    "position_1based": position_1based,
                    "off_target_sequence": off_target_sequence,
                    "strand": strand,
                    "mismatch_count": mismatch_count,
                    "raw_line": line,
                })

            except Exception as e:
                rows.append({
                    "accession": accession,
                    "source_file": str(path),
                    "line_num": line_num,
                    "parse_status": f"FAILED_{type(e).__name__}",
                    "raw_line": line,
                })

    return rows

def accession_from_filename(path: Path):
    return path.name.replace("_expanded_offtargets.txt", "")

def main():
    guides = load_guides()
    manifest = pd.read_csv(MANIFEST_FILE)

    raw_files = sorted(RAW_DIR.glob("*_expanded_offtargets.txt"))

    all_rows = []
    for raw_file in raw_files:
        accession = accession_from_filename(raw_file)
        all_rows.extend(parse_raw_file(raw_file, accession))

    hits = pd.DataFrame(all_rows)

    if hits.empty:
        raise ValueError("No rows parsed from raw Cas-OFFinder outputs.")

    passed = hits[hits["parse_status"] == "PASS"].copy()
    failed = hits[hits["parse_status"] != "PASS"].copy()

    passed = passed.merge(
        guides,
        left_on="query_sequence",
        right_on="casoffinder_query",
        how="left"
    )

    passed["guide_mapping_status"] = passed["guide_id"].apply(
        lambda x: "MAPPED" if pd.notna(x) else "UNMAPPED"
    )

    passed = passed.merge(
        manifest[["accession", "genome_file", "config_file", "raw_output_file"]],
        on="accession",
        how="left"
    )

    passed.to_csv(OUT_HITS, index=False)

    summary = []
    summary.append("Expanded Cas-OFFinder Parsed Summary")
    summary.append("===================================")
    summary.append(f"Raw files found: {len(raw_files)}")
    summary.append(f"Expected raw files from manifest: {manifest.shape[0]}")
    summary.append(f"Total raw lines parsed: {hits.shape[0]}")
    summary.append(f"PASS parsed rows: {passed.shape[0]}")
    summary.append(f"Failed parsed rows: {failed.shape[0]}")
    summary.append(f"Unique accessions with hits: {passed['accession'].nunique()}")
    summary.append(f"Unique mapped guides: {passed['guide_id'].nunique(dropna=True)}")
    summary.append(f"Unmapped rows: {(passed['guide_mapping_status'] == 'UNMAPPED').sum()}")
    summary.append("")
    summary.append("Hit count by target gene:")
    summary.append(str(passed["target_gene"].value_counts(dropna=False)))
    summary.append("")
    summary.append("Hit count by mismatch count:")
    summary.append(str(passed["mismatch_count"].value_counts(dropna=False).sort_index()))
    summary.append("")
    summary.append("Top 10 guides by expanded hit count:")
    summary.append(str(passed["guide_id"].value_counts(dropna=False).head(10)))
    summary.append("")
    summary.append(f"Output hit table: {OUT_HITS}")
    summary.append("QC status: " + ("PASS" if len(raw_files) == manifest.shape[0] and failed.shape[0] == 0 and (passed['guide_mapping_status'] == 'UNMAPPED').sum() == 0 else "CHECK"))

    OUT_SUMMARY.write_text("\n".join(summary) + "\n")
    print("\n".join(summary))

if __name__ == "__main__":
    main()
