from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "results_external"
    / "crisot"
    / "merged"
    / "crisot_master_merged.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "results_external"
    / "crisot"
    / "merged"
    / "crisot_guide_level_summary.csv"
)


def main() -> None:
    df = pd.read_csv(INPUT_FILE)

    required_columns = {
        "spacer",
        "pam",
        "guide_23mer",
        "Mutation",
        "CRISOT-Score",
        "CRISOT-Spec",
        "delta_Spec",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    wt = df[df["Mutation"].astype(str).str.upper() == "WT"].copy()
    mutants = df[df["Mutation"].astype(str).str.upper() != "WT"].copy()

    wt_summary = wt[
        [
            "spacer",
            "pam",
            "guide_23mer",
            "CRISOT-Score",
            "CRISOT-Spec",
        ]
    ].rename(
        columns={
            "CRISOT-Score": "crisot_wt_score",
            "CRISOT-Spec": "crisot_wt_specificity",
        }
    )

    mutant_summary = (
        mutants.groupby(["spacer", "pam", "guide_23mer"])
        .agg(
            crisot_mutation_count=("Mutation", "count"),
            crisot_mean_mutant_score=("CRISOT-Score", "mean"),
            crisot_min_mutant_score=("CRISOT-Score", "min"),
            crisot_max_mutant_score=("CRISOT-Score", "max"),
            crisot_mean_delta_spec=("delta_Spec", "mean"),
            crisot_max_delta_spec=("delta_Spec", "max"),
        )
        .reset_index()
    )

    summary = wt_summary.merge(
        mutant_summary,
        on=["spacer", "pam", "guide_23mer"],
        how="left",
    )

    numeric_cols = [
        "crisot_wt_score",
        "crisot_wt_specificity",
        "crisot_mean_mutant_score",
        "crisot_min_mutant_score",
        "crisot_max_mutant_score",
        "crisot_mean_delta_spec",
        "crisot_max_delta_spec",
    ]

    for col in numeric_cols:
        summary[col] = summary[col].round(6)

    summary.to_csv(OUTPUT_FILE, index=False)

    print(f"Input file: {INPUT_FILE}")
    print(f"Output file: {OUTPUT_FILE}")
    print(f"Guide-level CRISOT rows: {len(summary)}")
    print()
    print(summary.head(10).to_string(index=False))


if __name__ == "__main__":
    main()