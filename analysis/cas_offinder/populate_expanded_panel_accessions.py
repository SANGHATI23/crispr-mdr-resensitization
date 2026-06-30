#!/usr/bin/env python3

import json
import subprocess
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(".").resolve()

OUT_PATH = PROJECT_ROOT / "data_cas_offinder/expanded_panel/metadata/expanded_panel_accession_candidates.csv"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

COLUMNS = [
    "target_gene",
    "panel_category",
    "organism_or_element",
    "accession",
    "assembly_name",
    "assembly_level",
    "source_database",
    "strain_name",
    "target_gene_confirmed",
    "plasmid_or_element_type",
    "include_in_final_panel",
    "reason_for_inclusion",
    "notes",
]

SEARCH_PLAN = [
    {
        "target_gene": "mecA",
        "panel_category": "complete_chromosome",
        "organism_or_element": "Staphylococcus aureus MRSA",
        "taxon": "Staphylococcus aureus",
        "limit": 15,
        "gene_hint": "mecA",
        "plasmid_or_element_type": "SCCmec/MRSA chromosome",
    },
    {
        "target_gene": "blaKPC",
        "panel_category": "complete_chromosome",
        "organism_or_element": "Klebsiella pneumoniae",
        "taxon": "Klebsiella pneumoniae",
        "limit": 10,
        "gene_hint": "KPC",
        "plasmid_or_element_type": "chromosome",
    },
    {
        "target_gene": "blaNDM1",
        "panel_category": "complete_chromosome",
        "organism_or_element": "Klebsiella pneumoniae",
        "taxon": "Klebsiella pneumoniae",
        "limit": 8,
        "gene_hint": "NDM",
        "plasmid_or_element_type": "chromosome",
    },
    {
        "target_gene": "blaNDM1",
        "panel_category": "complete_chromosome",
        "organism_or_element": "Escherichia coli",
        "taxon": "Escherichia coli",
        "limit": 8,
        "gene_hint": "NDM",
        "plasmid_or_element_type": "chromosome",
    },
    {
        "target_gene": "mcr1",
        "panel_category": "complete_chromosome",
        "organism_or_element": "Escherichia coli",
        "taxon": "Escherichia coli",
        "limit": 10,
        "gene_hint": "mcr",
        "plasmid_or_element_type": "chromosome",
    },
]


def run_cmd(cmd):
    print("\nRunning:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        raise RuntimeError("Command failed")
    return result.stdout


def safe_get(d, path, default=""):
    cur = d
    for key in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key)
    return cur if cur is not None else default


def search_assemblies(taxon, limit):
    cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        taxon,
        "--assembly-level",
        "complete",
        "--limit",
        str(limit),
    ]
    stdout = run_cmd(cmd)
    data = json.loads(stdout)
    return data.get("reports", [])


def normalize_report(report, plan):
    accession = safe_get(report, ["accession"])
    assembly_name = safe_get(report, ["assembly_info", "assembly_name"])
    assembly_level = safe_get(report, ["assembly_info", "assembly_level"])
    source_database = safe_get(report, ["source_database"])
    strain_name = safe_get(report, ["organism", "infraspecific_names", "strain"])

    return {
        "target_gene": plan["target_gene"],
        "panel_category": plan["panel_category"],
        "organism_or_element": plan["organism_or_element"],
        "accession": accession,
        "assembly_name": assembly_name,
        "assembly_level": assembly_level,
        "source_database": source_database or "NCBI",
        "strain_name": strain_name,
        "target_gene_confirmed": "Needs confirmation",
        "plasmid_or_element_type": plan["plasmid_or_element_type"],
        "include_in_final_panel": "Candidate",
        "reason_for_inclusion": f"Candidate complete genome from NCBI Datasets search; target hint={plan['gene_hint']}",
        "notes": "Target gene presence must be confirmed by AMRFinder/BLAST/annotation before final claim",
    }


def main():
    existing = pd.DataFrame(columns=COLUMNS)
    if OUT_PATH.exists():
        existing = pd.read_csv(OUT_PATH)
        for col in COLUMNS:
            if col not in existing.columns:
                existing[col] = ""
        existing = existing[COLUMNS]

    new_rows = []

    for plan in SEARCH_PLAN:
        reports = search_assemblies(plan["taxon"], plan["limit"])
        print(f"Found {len(reports)} reports for {plan['target_gene']} / {plan['taxon']}")

        for report in reports:
            row = normalize_report(report, plan)
            if row["accession"]:
                new_rows.append(row)

    new_df = pd.DataFrame(new_rows, columns=COLUMNS)

    combined = pd.concat([existing, new_df], ignore_index=True)
    combined = combined.drop_duplicates(subset=["target_gene", "panel_category", "accession"], keep="first")
    combined = combined[COLUMNS]

    combined.to_csv(OUT_PATH, index=False)

    print("\nWROTE:", OUT_PATH)
    print("Rows:", combined.shape[0])
    print("Columns:", combined.shape[1])
    print("\nTarget gene counts:")
    print(combined["target_gene"].value_counts(dropna=False))
    print("\nPanel category counts:")
    print(combined["panel_category"].value_counts(dropna=False))
    print("\nInclude in final panel counts:")
    print(combined["include_in_final_panel"].value_counts(dropna=False))


if __name__ == "__main__":
    main()