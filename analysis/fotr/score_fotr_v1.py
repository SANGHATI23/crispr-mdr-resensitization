from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results_rna" / "all_guides_rna_structure.csv"
OUTPUT_DIR = PROJECT_ROOT / "results_fotr"
OUTPUT_FILE = OUTPUT_DIR / "all_guides_fotr_scores.csv"


def minmax_normalize(series: pd.Series) -> pd.Series:
    min_value = series.min()
    max_value = series.max()

    if max_value == min_value:
        return pd.Series([0.0] * len(series), index=series.index)

    return (series - min_value) / (max_value - min_value)


def compute_binding_risk(df: pd.DataFrame) -> pd.Series:
    """
    Higher value = higher predicted binding/off-target risk.
    This combines mismatch-tier hit burden and inverse specificity.
    """
    weighted_hits = (
        df["offtarget_hits_0mm"] * 1.00
        + df["offtarget_hits_1mm"] * 0.75
        + df["offtarget_hits_2mm"] * 0.50
        + df["offtarget_hits_3mm"] * 0.25
    )

    hit_risk = minmax_normalize(weighted_hits)
    inverse_specificity_risk = 1.0 - (df["specificity_score"] / 100.0)

    binding_risk = (0.60 * hit_risk) + (0.40 * inverse_specificity_risk)

    return binding_risk.clip(0, 1)


def compute_rna_structure_risk(df: pd.DataFrame) -> pd.Series:
    """
    Higher value = higher predicted RNA structural occlusion risk.
    Seed pairing is weighted more strongly than total spacer pairing.
    """
    rna_risk = (
        0.60 * df["seed_paired_fraction"]
        + 0.40 * df["spacer_paired_fraction"]
    )

    return rna_risk.clip(0, 1)


def main() -> None:
    df = pd.read_csv(INPUT_FILE)

    required_columns = {
        "gene",
        "spacer",
        "pam",
        "final_score",
        "specificity_score",
        "conservation_score",
        "offtarget_hits_0mm",
        "offtarget_hits_1mm",
        "offtarget_hits_2mm",
        "offtarget_hits_3mm",
        "seed_paired_fraction",
        "spacer_paired_fraction",
        "spacer_accessibility_score",
        "rna_accessibility_class",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    # Existing guide utility features
    df["activity_component"] = df["final_score"] / 100.0
    df["specificity_component"] = df["specificity_score"] / 100.0
    df["conservation_component"] = df["conservation_score"] / 100.0
    df["rna_accessibility_component"] = df["spacer_accessibility_score"]

    # Risk features
    df["binding_risk"] = compute_binding_risk(df)
    df["rna_structure_risk"] = compute_rna_structure_risk(df)

    # Neutral placeholders for annotation-based modules.
    # These will be upgraded later using AMR/virulence/essential-gene/context annotation.
    df["functional_severity"] = 1.0
    df["context_penalty"] = 1.0

    df["fotr_risk_score"] = (
        df["binding_risk"]
        * df["functional_severity"]
        * df["conservation_component"]
        * df["context_penalty"]
    )

    # Final prioritization score:
    # Higher = better guide candidate after activity, specificity, conservation, RNA accessibility, and FOTR risk.
    df["fotr_priority_score"] = (
        0.30 * df["activity_component"]
        + 0.25 * df["specificity_component"]
        + 0.25 * df["conservation_component"]
        + 0.20 * df["rna_accessibility_component"]
        - 0.20 * df["fotr_risk_score"]
        - 0.15 * df["rna_structure_risk"]
    )

    df["fotr_priority_score"] = (df["fotr_priority_score"] * 100).round(3)
    df["fotr_risk_score"] = (df["fotr_risk_score"] * 100).round(3)
    df["binding_risk"] = (df["binding_risk"] * 100).round(3)
    df["rna_structure_risk"] = (df["rna_structure_risk"] * 100).round(3)

    df = df.sort_values(
        by=["gene", "fotr_priority_score"],
        ascending=[True, False],
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Input file: {INPUT_FILE}")
    print(f"Output file: {OUTPUT_FILE}")
    print(f"Total guides scored: {len(df)}")
    print()
    print("Top 3 FOTR-prioritized guides per gene:")

    cols = [
        "gene",
        "position",
        "strand",
        "spacer",
        "pam",
        "final_score",
        "specificity_score",
        "conservation_score",
        "rna_accessibility_class",
        "binding_risk",
        "rna_structure_risk",
        "fotr_risk_score",
        "fotr_priority_score",
    ]

    print(df.groupby("gene").head(3)[cols].to_string(index=False))


if __name__ == "__main__":
    main()