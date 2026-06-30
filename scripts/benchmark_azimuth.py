from pathlib import Path
import sys
import pandas as pd
import numpy as np

try:
    import azimuth.model_comparison as azimuth_model
except Exception as e:
    raise RuntimeError(
        "\nFailed to import azimuth.model_comparison.\n"
        "This usually means the installed 'azimuth' package is the old Python-2-style version.\n"
        "Use a Python-3-compatible Azimuth installation, then rerun this script.\n"
        f"\nOriginal error:\n{e}\n"
    )


REPO_ROOT = Path(__file__).resolve().parents[1]

INPUT_CSV = REPO_ROOT / "results_panstrain" / "all_panstrain_guide_candidates.csv"
OUTPUT_DIR = REPO_ROOT / "results_external" / "azimuth"
OUTPUT_CSV = OUTPUT_DIR / "all_panstrain_guides_with_azimuth.csv"
TOP20_CSV = OUTPUT_DIR / "top20_azimuth_guides.csv"
SUMMARY_CSV = OUTPUT_DIR / "azimuth_summary.csv"


def build_30mer(spacer: str, pam: str) -> str:
    """
    Build a simple 30mer for Azimuth:
    4 nt left padding + 20 nt spacer + 3 nt PAM + 3 nt right padding = 30 nt

    This is a placeholder context for first-pass benchmarking.
    """
    spacer = str(spacer).strip().upper()
    pam = str(pam).strip().upper()

    left_pad = "NNNN"
    right_pad = "NNN"

    seq = f"{left_pad}{spacer}{pam}{right_pad}"

    if len(spacer) != 20:
        raise ValueError(f"Spacer must be length 20, got {len(spacer)}: {spacer}")

    if len(pam) != 3:
        raise ValueError(f"PAM must be length 3, got {len(pam)}: {pam}")

    if len(seq) != 30:
        raise ValueError(f"Invalid 30mer length {len(seq)} for sequence {seq}")

    return seq


def run_azimuth(seqs: list[str]) -> np.ndarray:
    """
    Run Azimuth in no-position mode.
    """
    seq_array = np.array(seqs, dtype=object)

    preds = azimuth_model.predict(
        seq_array,
        None,   # aa_cut
        None    # percent_peptide
    )

    preds = np.array(preds, dtype=float)
    preds = np.clip(preds, 0.0, 1.0)
    return preds


def validate_input_columns(df: pd.DataFrame) -> None:
    required_cols = ["gene", "position", "spacer", "pam"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def main() -> None:
    print("\n=== Azimuth Benchmark ===")

    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    validate_input_columns(df)

    df = df.copy()

    print(f"Loaded input rows: {len(df)}")
    print(f"Input file: {INPUT_CSV}")

    df["azimuth_30mer"] = df.apply(
        lambda row: build_30mer(row["spacer"], row["pam"]),
        axis=1
    )

    print("Built Azimuth 30mer sequences.")

    df["azimuth_raw"] = run_azimuth(df["azimuth_30mer"].tolist())
    df["azimuth_score_100"] = (df["azimuth_raw"] * 100).round(2)

    df["rank_azimuth"] = (
        df["azimuth_score_100"]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df.to_csv(OUTPUT_CSV, index=False)

    sort_cols = ["azimuth_score_100"]
    ascending = [False]

    if "specificity_score" in df.columns:
        sort_cols.append("specificity_score")
        ascending.append(False)

    if "conservation_score" in df.columns:
        sort_cols.append("conservation_score")
        ascending.append(False)

    top20 = df.sort_values(sort_cols, ascending=ascending).head(20)
    top20.to_csv(TOP20_CSV, index=False)

    summary = pd.DataFrame([
        {
            "n_guides": len(df),
            "azimuth_mean": round(df["azimuth_score_100"].mean(), 4),
            "azimuth_median": round(df["azimuth_score_100"].median(), 4),
            "azimuth_min": round(df["azimuth_score_100"].min(), 4),
            "azimuth_max": round(df["azimuth_score_100"].max(), 4),
        }
    ])
    summary.to_csv(SUMMARY_CSV, index=False)

    print("\nAzimuth benchmark completed successfully.")
    print(f"Saved full output: {OUTPUT_CSV}")
    print(f"Saved top 20:      {TOP20_CSV}")
    print(f"Saved summary:     {SUMMARY_CSV}")

    cols_to_show = ["gene", "position", "spacer", "pam", "azimuth_score_100", "rank_azimuth"]
    cols_to_show = [c for c in cols_to_show if c in top20.columns]

    print("\nTop 10 Azimuth-ranked guides:")
    print(top20[cols_to_show].head(10).to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERROR: {exc}")
        sys.exit(1)