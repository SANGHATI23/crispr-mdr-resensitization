#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def minmax(s):
    x = pd.to_numeric(s, errors="coerce")
    if x.max() == x.min():
        return pd.Series(0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", default="results_external")
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    required = ["final_score", "conservation", "CFD_score"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df["norm_final_score"] = minmax(df["final_score"])
    df["norm_conservation"] = minmax(df["conservation"])
    df["norm_CFD_score"] = minmax(df["CFD_score"])

    df["external_strain_aware_score"] = (
        0.70 * df["norm_final_score"] +
        0.20 * df["norm_conservation"] +
        0.10 * df["norm_CFD_score"]
    )

    df["internal_rank"] = df["final_score"].rank(ascending=False, method="min").astype(int)
    df["external_strain_aware_rank"] = df["external_strain_aware_score"].rank(
        ascending=False, method="min"
    ).astype(int)

    df["rank_shift_external_minus_internal"] = (
        df["external_strain_aware_rank"] - df["internal_rank"]
    )

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    df.sort_values("external_strain_aware_rank").to_csv(
        outdir / "strain_aware_comparison.csv", index=False
    )

    rows = []
    for k in [10, 20, 50]:
        top_internal = set(df.sort_values("internal_rank").head(k).index)
        top_external = set(df.sort_values("external_strain_aware_rank").head(k).index)
        overlap = len(top_internal & top_external)

        rows.append({
            "top_k": k,
            "overlap": overlap,
            "percent": round(overlap / k * 100, 2)
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "strain_overlap_summary.csv", index=False)

    print("\nStrain-aware baseline complete\n")
    print(summary.to_string(index=False))

    print("\nInterpretation:")
    print("External strain-aware ranking = 70% internal composite + 20% conservation + 10% CFD.")
    print("This tests whether guide prioritization remains stable after conservation is explicitly upweighted.")


if __name__ == "__main__":
    main()
