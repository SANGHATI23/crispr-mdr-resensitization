#!/usr/bin/env python3
"""
Create the planned expanded Cas-OFFinder genome-panel manifest.

The manifest defines the genome and mobile-element categories needed to
complete the genome-wide off-target analysis for:

    blaKPC
    blaNDM1
    mcr1
    mecA

This script does not download genomes. It creates an auditable design
document specifying the minimum panel composition required before the
expanded Cas-OFFinder run.

Outputs
-------
data_cas_offinder/expanded_panel/metadata/
    expanded_genome_panel_manifest.csv
    expanded_genome_panel_requirements_summary.txt
"""

from pathlib import Path
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data_cas_offinder"
    / "expanded_panel"
    / "metadata"
)

MANIFEST_FILE = OUTPUT_DIR / "expanded_genome_panel_manifest.csv"

SUMMARY_FILE = (
    OUTPUT_DIR
    / "expanded_genome_panel_requirements_summary.txt"
)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = [
        # -----------------------------------------------------
        # blaKPC panel
        # -----------------------------------------------------
        {
            "target_gene": "blaKPC",
            "organism_or_element": "Klebsiella pneumoniae",
            "panel_category": "complete_chromosome",
            "required_count": 10,
            "priority": "Required",
            "selection_requirement": (
                "Clinically relevant KPC-positive complete or "
                "chromosome-level assemblies"
            ),
            "purpose": (
                "Measure chromosomal off-target burden across "
                "diverse K. pneumoniae backgrounds"
            ),
        },
        {
            "target_gene": "blaKPC",
            "organism_or_element": "Enterobacterales",
            "panel_category": "cross_species_genome",
            "required_count": 5,
            "priority": "Required",
            "selection_requirement": (
                "KPC-positive Enterobacter, Escherichia, Citrobacter, "
                "or Serratia assemblies"
            ),
            "purpose": (
                "Test whether off-target sites occur outside "
                "K. pneumoniae"
            ),
        },
        {
            "target_gene": "blaKPC",
            "organism_or_element": "KPC plasmids",
            "panel_category": "complete_plasmid",
            "required_count": 10,
            "priority": "Required",
            "selection_requirement": (
                "Complete blaKPC-bearing plasmids representing common "
                "replicon backgrounds"
            ),
            "purpose": (
                "Capture plasmid-specific and mobile-element "
                "off-target sites"
            ),
        },

        # -----------------------------------------------------
        # blaNDM1 panel
        # -----------------------------------------------------
        {
            "target_gene": "blaNDM1",
            "organism_or_element": "Klebsiella pneumoniae",
            "panel_category": "complete_chromosome",
            "required_count": 8,
            "priority": "Required",
            "selection_requirement": (
                "NDM-positive clinically relevant complete or "
                "chromosome-level assemblies"
            ),
            "purpose": (
                "Measure NDM-guide off-target burden in Klebsiella "
                "backgrounds"
            ),
        },
        {
            "target_gene": "blaNDM1",
            "organism_or_element": "Escherichia coli",
            "panel_category": "complete_chromosome",
            "required_count": 8,
            "priority": "Required",
            "selection_requirement": (
                "NDM-positive complete or chromosome-level "
                "E. coli assemblies"
            ),
            "purpose": (
                "Measure NDM-guide off-target burden in E. coli "
                "backgrounds"
            ),
        },
        {
            "target_gene": "blaNDM1",
            "organism_or_element": "NDM plasmids and integrons",
            "panel_category": "complete_plasmid_or_integron",
            "required_count": 10,
            "priority": "Required",
            "selection_requirement": (
                "Complete blaNDM1-bearing plasmids or integron-containing "
                "sequences from diverse replicon types"
            ),
            "purpose": (
                "Capture plasmid, integron, and mobile-element "
                "off-target contexts"
            ),
        },

        # -----------------------------------------------------
        # mcr1 panel
        # -----------------------------------------------------
        {
            "target_gene": "mcr1",
            "organism_or_element": "Escherichia coli",
            "panel_category": "complete_chromosome",
            "required_count": 10,
            "priority": "Required",
            "selection_requirement": (
                "mcr-1-positive complete or chromosome-level "
                "E. coli assemblies"
            ),
            "purpose": (
                "Measure off-target burden across representative "
                "mcr-1 host backgrounds"
            ),
        },
        {
            "target_gene": "mcr1",
            "organism_or_element": "Other Enterobacterales",
            "panel_category": "cross_species_genome",
            "required_count": 5,
            "priority": "Recommended",
            "selection_requirement": (
                "mcr-1-positive Klebsiella, Salmonella, or other "
                "Enterobacterales assemblies"
            ),
            "purpose": (
                "Evaluate cross-species conservation of off-target risk"
            ),
        },
        {
            "target_gene": "mcr1",
            "organism_or_element": "mcr-1 plasmids",
            "panel_category": "complete_plasmid",
            "required_count": 12,
            "priority": "Required",
            "selection_requirement": (
                "Complete mcr-1 plasmids including IncI2, IncX4, "
                "IncHI2, and other major replicon contexts"
            ),
            "purpose": (
                "Measure plasmid-context off-target burden and "
                "mobile-element recurrence"
            ),
        },

        # -----------------------------------------------------
        # mecA panel
        # -----------------------------------------------------
        {
            "target_gene": "mecA",
            "organism_or_element": "Staphylococcus aureus MRSA",
            "panel_category": "complete_chromosome",
            "required_count": 15,
            "priority": "Required",
            "selection_requirement": (
                "mecA-positive complete MRSA genomes representing "
                "diverse clonal and SCCmec backgrounds"
            ),
            "purpose": (
                "Measure mecA-guide off-target burden across "
                "clinically relevant MRSA genomes"
            ),
        },
        {
            "target_gene": "mecA",
            "organism_or_element": "Coagulase-negative Staphylococci",
            "panel_category": "cross_species_genome",
            "required_count": 5,
            "priority": "Recommended",
            "selection_requirement": (
                "mecA-positive complete Staphylococcus epidermidis "
                "or related species assemblies"
            ),
            "purpose": (
                "Assess cross-species risk in mecA reservoirs"
            ),
        },
        {
            "target_gene": "mecA",
            "organism_or_element": "SCCmec elements",
            "panel_category": "complete_mobile_element",
            "required_count": 10,
            "priority": "Required",
            "selection_requirement": (
                "Representative SCCmec types containing mecA and "
                "associated recombinase/resistance regions"
            ),
            "purpose": (
                "Map guide matches and off-targets directly within "
                "SCCmec contexts"
            ),
        },

        # -----------------------------------------------------
        # Negative/reference controls
        # -----------------------------------------------------
        {
            "target_gene": "all",
            "organism_or_element": "Target-negative reference genomes",
            "panel_category": "negative_control_genome",
            "required_count": 8,
            "priority": "Required",
            "selection_requirement": (
                "High-quality genomes without the corresponding "
                "target resistance determinant"
            ),
            "purpose": (
                "Estimate background off-target burden unrelated "
                "to target-gene presence"
            ),
        },
    ]

    manifest = pd.DataFrame(rows)

    manifest.insert(
        0,
        "panel_record_id",
        [
            f"panel_{index:02d}"
            for index in range(1, len(manifest) + 1)
        ],
    )

    manifest["download_status"] = "Not_started"
    manifest["accession"] = ""
    manifest["assembly_level"] = ""
    manifest["strain_name"] = ""
    manifest["sequence_type_or_lineage"] = ""
    manifest["plasmid_replicon_or_element_type"] = ""
    manifest["source_database"] = ""
    manifest["fasta_filename"] = ""
    manifest["gff_filename"] = ""
    manifest["target_gene_confirmed"] = ""
    manifest["quality_check_status"] = ""
    manifest["notes"] = ""

    manifest.to_csv(MANIFEST_FILE, index=False)

    required_total = int(
        manifest.loc[
            manifest["priority"] == "Required",
            "required_count",
        ].sum()
    )

    recommended_total = int(
        manifest.loc[
            manifest["priority"] == "Recommended",
            "required_count",
        ].sum()
    )

    total_planned = int(manifest["required_count"].sum())

    summary_lines = [
        "Expanded Cas-OFFinder Genome Panel Requirements",
        "=" * 72,
        "",
        f"Manifest records: {len(manifest)}",
        f"Required sequences/genomes: {required_total}",
        f"Recommended additional sequences/genomes: {recommended_total}",
        f"Total planned panel size: {total_planned}",
        "",
        "Planned count by target gene",
        "-" * 72,
        manifest.groupby("target_gene")[
            "required_count"
        ].sum().to_string(),
        "",
        "Planned count by panel category",
        "-" * 72,
        manifest.groupby("panel_category")[
            "required_count"
        ].sum().to_string(),
        "",
        "Important design rule",
        "-" * 72,
        (
            "The final panel must contain both bacterial chromosome/"
            "assembly sequences and complete plasmid or mobile-element "
            "sequences. A chromosome-only panel is insufficient for "
            "blaKPC, blaNDM1, mcr1, and mecA because these resistance "
            "determinants frequently occur in mobile genetic contexts."
        ),
        "",
        "Next action",
        "-" * 72,
        (
            "Populate accession, strain, assembly quality, target-gene "
            "confirmation, plasmid replicon, SCCmec type, FASTA filename, "
            "and GFF filename before running the expanded Cas-OFFinder "
            "analysis."
        ),
    ]

    SUMMARY_FILE.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print("=" * 72)
    print("Expanded genome-panel manifest created")
    print("=" * 72)
    print(f"Manifest rows: {len(manifest)}")
    print(f"Required panel count: {required_total}")
    print(f"Recommended additional count: {recommended_total}")
    print(f"Total planned panel count: {total_planned}")
    print("\nFiles written:")
    print(f"- {MANIFEST_FILE.relative_to(PROJECT_ROOT)}")
    print(f"- {SUMMARY_FILE.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()