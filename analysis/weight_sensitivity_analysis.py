import os
import pandas as pd

# =========================
# INPUT FILES
# =========================
PHASE1_INPUT = "results/all_guide_candidates.csv"
PHASE2_INPUT = "results_panstrain/all_panstrain_guide_candidates.csv"

OUTDIR = "weight_sensitivity_results"
os.makedirs(OUTDIR, exist_ok=True)


# =========================
# CLASSIFICATION FUNCTIONS
# =========================
def classify_phase1(score):
    if score >= 85:
        return "Excellent"
    elif score >= 70:
        return "Good"
    elif score >= 50:
        return "Moderate"
    else:
        return "Poor"


def classify_phase2(score):
    if score >= 85:
        return "Excellent"
    elif score >= 70:
        return "Good"
    elif score >= 50:
        return "Moderate"
    else:
        return "Poor"


# =========================
# HELPER FUNCTIONS
# =========================
def jaccard_top_n(df_a, df_b, n=20):
    a = set(df_a.sort_values("final_score_new", ascending=False).head(n)["guide_key"])
    b = set(df_b.sort_values("final_score_new", ascending=False).head(n)["guide_key"])
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return round(len(a & b) / len(a | b), 4)


def overlap_count_top_n(df_a, df_b, n=20):
    a = set(df_a.sort_values("final_score_new", ascending=False).head(n)["guide_key"])
    b = set(df_b.sort_values("final_score_new", ascending=False).head(n)["guide_key"])
    return len(a & b)


# =========================
# PHASE 1: ALL VARIATIONS
# On-target / Specificity
# =========================
def run_phase1_all_variations():
    df = pd.read_csv(PHASE1_INPUT)

    required = {"gene", "position", "spacer", "on_target_score", "specificity_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Phase 1 file missing columns: {missing}")

    df["guide_key"] = (
        df["gene"].astype(str) + "|" +
        df["position"].astype(str) + "|" +
        df["spacer"].astype(str)
    )

    # All 11 variations from 100/0 to 0/100
    weights = []
    for w_on_int in range(100, -1, -10):
        w_spec_int = 100 - w_on_int
        weights.append((w_on_int / 100.0, w_spec_int / 100.0))

    baseline = None
    summary_rows = []
    all_rows = []

    for w_on, w_spec in weights:
        tmp = df.copy()
        scheme = f"{int(w_on*100)}/{int(w_spec*100)}"
        tmp["weight_scheme"] = scheme

        tmp["final_score_new"] = (
            w_on * tmp["on_target_score"] +
            w_spec * tmp["specificity_score"]
        ).round(1)

        tmp["classification_new"] = tmp["final_score_new"].apply(classify_phase1)

        # Use 60/40 as baseline for comparison
        if (w_on, w_spec) == (0.60, 0.40):
            baseline = tmp.copy()

        all_rows.append(tmp)

    all_df = pd.concat(all_rows, ignore_index=True)

    if baseline is None:
        raise ValueError("60/40 baseline not found in Phase 1 weights.")

    for scheme in all_df["weight_scheme"].unique():
        tmp = all_df[all_df["weight_scheme"] == scheme].copy()
        top20 = tmp.sort_values("final_score_new", ascending=False).head(20)

        summary_rows.append({
            "weight_scheme": scheme,
            "mean_final_score": round(tmp["final_score_new"].mean(), 2),
            "median_final_score": round(tmp["final_score_new"].median(), 2),
            "best_final_score": round(tmp["final_score_new"].max(), 2),
            "excellent_count": int((tmp["classification_new"] == "Excellent").sum()),
            "good_count": int((tmp["classification_new"] == "Good").sum()),
            "moderate_count": int((tmp["classification_new"] == "Moderate").sum()),
            "poor_count": int((tmp["classification_new"] == "Poor").sum()),
            "top20_overlap_count_vs_60_40": overlap_count_top_n(baseline, tmp, n=20),
            "top20_jaccard_vs_60_40": jaccard_top_n(baseline, tmp, n=20),
            "top_guide_gene": top20.iloc[0]["gene"] if len(top20) else None,
            "top_guide_spacer": top20.iloc[0]["spacer"] if len(top20) else None,
            "top_guide_score": round(top20.iloc[0]["final_score_new"], 2) if len(top20) else None
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(
        by="weight_scheme",
        key=lambda s: s.str.split("/").str[0].astype(int),
        ascending=False
    )

    summary_df.to_csv(os.path.join(OUTDIR, "phase1_weight_summary_all_variations.csv"), index=False)
    all_df.to_csv(os.path.join(OUTDIR, "phase1_all_weighted_guides_all_variations.csv"), index=False)

    print("Saved Phase 1 all-variation sensitivity results.")


# =========================
# PHASE 2: ALL VARIATIONS
# On-target / Specificity / Conservation
# =========================
def run_phase2_all_variations():
    df = pd.read_csv(PHASE2_INPUT)

    required = {
        "gene", "position", "spacer",
        "on_target_score", "specificity_score", "conservation_score"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Phase 2 file missing columns: {missing}")

    df["guide_key"] = (
        df["gene"].astype(str) + "|" +
        df["position"].astype(str) + "|" +
        df["spacer"].astype(str)
    )

    # All 66 combinations where:
    # w_on + w_spec + w_cons = 1.0
    # each weight moves in 0.1 steps
    weights = []
    for w_on_int in range(0, 101, 10):
        for w_spec_int in range(0, 101 - w_on_int, 10):
            w_cons_int = 100 - w_on_int - w_spec_int
            weights.append((w_on_int / 100.0, w_spec_int / 100.0, w_cons_int / 100.0))

    baseline = None
    summary_rows = []
    all_rows = []

    for w_on, w_spec, w_cons in weights:
        tmp = df.copy()
        scheme = f"{int(w_on*100)}/{int(w_spec*100)}/{int(w_cons*100)}"
        tmp["weight_scheme"] = scheme

        tmp["final_score_new"] = (
            w_on * tmp["on_target_score"] +
            w_spec * tmp["specificity_score"] +
            w_cons * tmp["conservation_score"]
        ).round(1)

        tmp["classification_new"] = tmp["final_score_new"].apply(classify_phase2)

        # Use manuscript baseline 40/30/30
        if (w_on, w_spec, w_cons) == (0.40, 0.30, 0.30):
            baseline = tmp.copy()

        all_rows.append(tmp)

    all_df = pd.concat(all_rows, ignore_index=True)

    if baseline is None:
        raise ValueError("40/30/30 baseline not found in Phase 2 weights.")

    for scheme in all_df["weight_scheme"].unique():
        tmp = all_df[all_df["weight_scheme"] == scheme].copy()
        top20 = tmp.sort_values("final_score_new", ascending=False).head(20)

        summary_rows.append({
            "weight_scheme": scheme,
            "mean_final_score": round(tmp["final_score_new"].mean(), 2),
            "median_final_score": round(tmp["final_score_new"].median(), 2),
            "best_final_score": round(tmp["final_score_new"].max(), 2),
            "excellent_count": int((tmp["classification_new"] == "Excellent").sum()),
            "good_count": int((tmp["classification_new"] == "Good").sum()),
            "moderate_count": int((tmp["classification_new"] == "Moderate").sum()),
            "poor_count": int((tmp["classification_new"] == "Poor").sum()),
            "top20_overlap_count_vs_40_30_30": overlap_count_top_n(baseline, tmp, n=20),
            "top20_jaccard_vs_40_30_30": jaccard_top_n(baseline, tmp, n=20),
            "top_guide_gene": top20.iloc[0]["gene"] if len(top20) else None,
            "top_guide_spacer": top20.iloc[0]["spacer"] if len(top20) else None,
            "top_guide_score": round(top20.iloc[0]["final_score_new"], 2) if len(top20) else None
        })

    summary_df = pd.DataFrame(summary_rows)

    def sort_key(series):
        return (
            series.str.split("/").str[0].astype(int) * 10000 +
            series.str.split("/").str[1].astype(int) * 100 +
            series.str.split("/").str[2].astype(int)
        )

    summary_df = summary_df.sort_values(by="weight_scheme", key=sort_key, ascending=False)

    summary_df.to_csv(os.path.join(OUTDIR, "phase2_weight_summary_all_variations.csv"), index=False)
    all_df.to_csv(os.path.join(OUTDIR, "phase2_all_weighted_guides_all_variations.csv"), index=False)

    print("Saved Phase 2 all-variation sensitivity results.")


# =========================
# OPTIONAL: PER-GENE SUMMARIES
# =========================
def run_phase1_per_gene_summary():
    all_df = pd.read_csv(os.path.join(OUTDIR, "phase1_all_weighted_guides_all_variations.csv"))

    summary = (
        all_df.groupby(["weight_scheme", "gene"])
        .agg(
            guide_count=("guide_key", "count"),
            mean_final_score=("final_score_new", "mean"),
            best_final_score=("final_score_new", "max"),
            excellent_count=("classification_new", lambda x: (x == "Excellent").sum()),
            good_count=("classification_new", lambda x: (x == "Good").sum()),
            moderate_count=("classification_new", lambda x: (x == "Moderate").sum()),
            poor_count=("classification_new", lambda x: (x == "Poor").sum()),
        )
        .reset_index()
    )

    summary["mean_final_score"] = summary["mean_final_score"].round(2)
    summary["best_final_score"] = summary["best_final_score"].round(2)

    summary.to_csv(os.path.join(OUTDIR, "phase1_per_gene_summary_all_variations.csv"), index=False)
    print("Saved Phase 1 per-gene summary.")


def run_phase2_per_gene_summary():
    all_df = pd.read_csv(os.path.join(OUTDIR, "phase2_all_weighted_guides_all_variations.csv"))

    summary = (
        all_df.groupby(["weight_scheme", "gene"])
        .agg(
            guide_count=("guide_key", "count"),
            mean_final_score=("final_score_new", "mean"),
            best_final_score=("final_score_new", "max"),
            mean_conservation=("conservation_score", "mean"),
            best_conservation=("conservation_score", "max"),
            excellent_count=("classification_new", lambda x: (x == "Excellent").sum()),
            good_count=("classification_new", lambda x: (x == "Good").sum()),
            moderate_count=("classification_new", lambda x: (x == "Moderate").sum()),
            poor_count=("classification_new", lambda x: (x == "Poor").sum()),
        )
        .reset_index()
    )

    summary["mean_final_score"] = summary["mean_final_score"].round(2)
    summary["best_final_score"] = summary["best_final_score"].round(2)
    summary["mean_conservation"] = summary["mean_conservation"].round(2)
    summary["best_conservation"] = summary["best_conservation"].round(2)

    summary.to_csv(os.path.join(OUTDIR, "phase2_per_gene_summary_all_variations.csv"), index=False)
    print("Saved Phase 2 per-gene summary.")


# =========================
# MAIN
# =========================
if __name__ == "__main__":
    run_phase1_all_variations()
    run_phase2_all_variations()
    run_phase1_per_gene_summary()
    run_phase2_per_gene_summary()
    print(f"\nDone. Check folder: {OUTDIR}")