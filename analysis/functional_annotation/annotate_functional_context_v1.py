"""
Functional/context annotation v1 for CRISPR-MDR FOTR.

Purpose:
- Add first-pass biological annotation columns to the FOTR+CRISOT guide table.
- Replace placeholder-only interpretation with explicit target-gene and genomic-context labels.
- This is intentionally conservative: it annotates target-side functional importance first,
  not true genome-wide off-target functional annotation yet.

Input:
    results_fotr/all_guides_fotr_scores_with_crisot.csv

Output:
    results_functional_annotation/all_guides_functional_context_v1.csv
    results_functional_annotation/functional_context_summary_by_gene.csv
"""

from pathlib import Path
import pandas as pd


INPUT_PATH = Path("results_fotr/all_guides_fotr_scores_with_crisot.csv")
OUTPUT_DIR = Path("results_functional_annotation")
OUTPUT_GUIDE_PATH = OUTPUT_DIR / "all_guides_functional_context_v1.csv"
OUTPUT_SUMMARY_PATH = OUTPUT_DIR / "functional_context_summary_by_gene.csv"


# First-pass AMR target annotation.
# These are target-side annotations, not experimental validation claims.
GENE_ANNOTATION = {
    "blaKPC": {
        "amr_class": "beta_lactamase_carbapenemase",
        "resistance_mechanism": "carbapenem_hydrolysis",
        "clinical_relevance": "high",
        "target_functional_importance": 1.00,
        "target_annotation_note": "KPC carbapenemase target; disruption is expected to reduce carbapenem resistance in carrying backgrounds.",
    },
    "blaNDM1": {
        "amr_class": "metallo_beta_lactamase",
        "resistance_mechanism": "carbapenem_hydrolysis",
        "clinical_relevance": "high",
        "target_functional_importance": 1.00,
        "target_annotation_note": "NDM-1 metallo-beta-lactamase target; disruption is expected to reduce carbapenem resistance in carrying backgrounds.",
    },
    "mcr1": {
        "amr_class": "phosphoethanolamine_transferase",
        "resistance_mechanism": "lipid_A_modification_colistin_resistance",
        "clinical_relevance": "high",
        "target_functional_importance": 0.95,
        "target_annotation_note": "mcr-1 colistin resistance target; disruption may reduce polymyxin/colistin resistance in carrying backgrounds.",
    },
    "mecA": {
        "amr_class": "penicillin_binding_protein",
        "resistance_mechanism": "altered_beta_lactam_target_PBP2a",
        "clinical_relevance": "high",
        "target_functional_importance": 1.00,
        "target_annotation_note": "mecA encodes PBP2a; disruption is expected to reduce methicillin/beta-lactam resistance in MRSA-like backgrounds.",
    },
}


def normalize_gene_name(gene: str) -> str:
    """Normalize gene names for consistent annotation lookup."""
    if pd.isna(gene):
        return ""
    g = str(gene).strip()
    aliases = {
        "blaNDM-1": "blaNDM1",
        "blaNDM_1": "blaNDM1",
        "mcr-1": "mcr1",
        "MCR1": "mcr1",
        "MECA": "mecA",
    }
    return aliases.get(g, g)


def assign_target_region(position: int) -> str:
    """
    First-pass coding-region position bin.

    This does not require protein-domain coordinates yet.
    It gives reviewers a transparent first target-context feature.
    """
    if pd.isna(position):
        return "unknown"

    pos = int(position)

    if pos < 100:
        return "early_coding_region"
    elif pos < 500:
        return "middle_coding_region"
    else:
        return "late_coding_region"


def assign_disruption_priority(position: int) -> float:
    """
    Conservative first-pass target disruption priority.

    Early coding targets are prioritized because frameshift/indel disruption near
    the N-terminus is more likely to disable the protein product.
    """
    region = assign_target_region(position)

    if region == "early_coding_region":
        return 1.00
    elif region == "middle_coding_region":
        return 0.85
    elif region == "late_coding_region":
        return 0.70
    return 0.50


def assign_context_label(gene: str) -> str:
    """
    First-pass genomic context label.

    This is a target-side context label based on known AMR gene biology.
    It does not yet replace full plasmid/contig/mobile-element annotation.
    """
    gene = normalize_gene_name(gene)

    if gene in {"blaKPC", "blaNDM1", "mcr1"}:
        return "often_plasmid_or_mobile_AMR_context"
    elif gene == "mecA":
        return "chromosomal_SCCmec_AMR_context"
    return "unknown_context"


def assign_context_penalty_v1(context_label: str) -> float:
    """
    First-pass context weight.

    Higher value means stronger AMR-relevant target context.
    This is not an off-target penalty yet; it is target-context importance.
    """
    if context_label == "often_plasmid_or_mobile_AMR_context":
        return 1.00
    elif context_label == "chromosomal_SCCmec_AMR_context":
        return 1.00
    return 0.75


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    required_cols = {"gene", "position", "spacer", "pam", "fotr_priority_score"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["gene_normalized"] = df["gene"].apply(normalize_gene_name)

    annotation_rows = []
    for _, row in df.iterrows():
        gene = row["gene_normalized"]
        ann = GENE_ANNOTATION.get(gene, {})

        target_region = assign_target_region(row["position"])
        target_disruption_priority = assign_disruption_priority(row["position"])
        target_context_label = assign_context_label(gene)
        target_context_weight_v1 = assign_context_penalty_v1(target_context_label)

        target_functional_importance = ann.get("target_functional_importance", 0.75)

        # Functional severity v1 is target-side functional importance x disruption location priority.
        # It intentionally does not claim genome-wide off-target severity yet.
        functional_severity_v1 = round(
            float(target_functional_importance) * float(target_disruption_priority),
            4,
        )

        annotation_rows.append(
            {
                "amr_class": ann.get("amr_class", "unknown"),
                "resistance_mechanism": ann.get("resistance_mechanism", "unknown"),
                "clinical_relevance": ann.get("clinical_relevance", "unknown"),
                "target_region_bin": target_region,
                "target_disruption_priority": target_disruption_priority,
                "target_functional_importance": target_functional_importance,
                "functional_severity_v1": functional_severity_v1,
                "target_context_label": target_context_label,
                "target_context_weight_v1": target_context_weight_v1,
                "target_annotation_note": ann.get("target_annotation_note", "No curated note available."),
                "annotation_scope": "target_side_first_pass_not_genomewide_offtarget",
            }
        )

    ann_df = pd.DataFrame(annotation_rows)
    out = pd.concat([df, ann_df], axis=1)

    out.to_csv(OUTPUT_GUIDE_PATH, index=False)

    summary = (
        out.groupby(["gene_normalized", "amr_class", "target_region_bin"])
        .agg(
            guide_count=("spacer", "count"),
            mean_fotr_priority_score=("fotr_priority_score", "mean"),
            mean_functional_severity_v1=("functional_severity_v1", "mean"),
            mean_target_disruption_priority=("target_disruption_priority", "mean"),
        )
        .reset_index()
        .sort_values(["gene_normalized", "target_region_bin"])
    )

    summary.to_csv(OUTPUT_SUMMARY_PATH, index=False)

    print(f"Wrote: {OUTPUT_GUIDE_PATH}")
    print(f"Wrote: {OUTPUT_SUMMARY_PATH}")
    print("\nGuide counts by gene and target region:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()