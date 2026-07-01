#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score, balanced_accuracy_score, f1_score

DATA_PATH = "results_mlcb/tables/mlcb_guide_dataset_with_local_context_embedding_ids.csv"
EMB_PATH = "results_mlcb/embeddings/dnabert2_local_context_embeddings.npy"

OUT_TABLE = "results_mlcb/tables/dnabert_calibration_metrics.csv"
OUT_BINS = "results_mlcb/tables/dnabert_calibration_curve_bins.csv"
OUT_FIG = "results_mlcb/figures/figure_dnabert_calibration_reliability_diagram.png"
OUT_TXT = "results_mlcb/tables/dnabert_calibration_summary.txt"

LABEL_CANDIDATES = [
    "weak_functional_severity_label",
    "label_high_functional_offtarget_risk",
    "label_high_priority_by_score_top_quartile",
]

def find_label_col(df):
    for col in LABEL_CANDIDATES:
        if col in df.columns:
            return col
    raise ValueError(f"No label column found. Tried {LABEL_CANDIDATES}")

def expected_calibration_error(y_true, y_prob, n_bins=10):
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)

    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    bin_rows = []

    for i in range(n_bins):
        left = bins[i]
        right = bins[i + 1]

        if i == n_bins - 1:
            mask = (y_prob >= left) & (y_prob <= right)
        else:
            mask = (y_prob >= left) & (y_prob < right)

        n_bin = int(mask.sum())
        if n_bin == 0:
            bin_rows.append({
                "bin": i + 1,
                "bin_left": left,
                "bin_right": right,
                "n": 0,
                "mean_confidence": np.nan,
                "fraction_positive": np.nan,
                "abs_gap": np.nan,
            })
            continue

        mean_conf = float(y_prob[mask].mean())
        frac_pos = float(y_true[mask].mean())
        gap = abs(frac_pos - mean_conf)

        ece += (n_bin / len(y_true)) * gap

        bin_rows.append({
            "bin": i + 1,
            "bin_left": left,
            "bin_right": right,
            "n": n_bin,
            "mean_confidence": mean_conf,
            "fraction_positive": frac_pos,
            "abs_gap": gap,
        })

    return float(ece), pd.DataFrame(bin_rows)

def evaluate_probs(y_true, y_prob, name):
    y_pred = (y_prob >= 0.5).astype(int)

    return {
        "model": name,
        "brier_score": brier_score_loss(y_true, y_prob),
        "log_loss": log_loss(y_true, y_prob, labels=[0, 1]),
        "ece_10_bins": expected_calibration_error(y_true, y_prob, n_bins=10)[0],
        "roc_auc": roc_auc_score(y_true, y_prob),
        "balanced_accuracy_at_0.5": balanced_accuracy_score(y_true, y_pred),
        "macro_f1_at_0.5": f1_score(y_true, y_pred, average="macro", zero_division=0),
    }

def main():
    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(f"Missing metadata table: {DATA_PATH}")
    if not os.path.exists(EMB_PATH):
        raise FileNotFoundError(f"Missing embedding matrix: {EMB_PATH}")

    df = pd.read_csv(DATA_PATH)
    X = np.load(EMB_PATH)

    label_col = find_label_col(df)
    y = df[label_col].astype(int).values

    print("Loaded metadata:", df.shape)
    print("Loaded embeddings:", X.shape)
    print("Label column:", label_col)
    print("Label counts:")
    print(pd.Series(y).value_counts())

    # We use cross_val_predict so every guide receives an out-of-fold probability.
    # This gives calibration diagnostics without evaluating on training predictions.
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    base_model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            solver="liblinear",
            random_state=42,
        )
    )

    platt_model = make_pipeline(
        StandardScaler(),
        CalibratedClassifierCV(
            estimator=LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42,
            ),
            method="sigmoid",
            cv=3,
        )
    )

    isotonic_model = make_pipeline(
        StandardScaler(),
        CalibratedClassifierCV(
            estimator=LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42,
            ),
            method="isotonic",
            cv=3,
        )
    )

    print("\nGenerating out-of-fold probabilities...")
    prob_uncal = cross_val_predict(base_model, X, y, cv=cv, method="predict_proba")[:, 1]
    prob_platt = cross_val_predict(platt_model, X, y, cv=cv, method="predict_proba")[:, 1]
    prob_iso = cross_val_predict(isotonic_model, X, y, cv=cv, method="predict_proba")[:, 1]

    metrics = pd.DataFrame([
        evaluate_probs(y, prob_uncal, "DNABERT logistic probe, uncalibrated"),
        evaluate_probs(y, prob_platt, "DNABERT logistic probe + Platt scaling"),
        evaluate_probs(y, prob_iso, "DNABERT logistic probe + isotonic regression"),
    ])

    metrics.to_csv(OUT_TABLE, index=False)

    # Save bin-level values for all methods
    bin_tables = []
    for name, prob in [
        ("uncalibrated", prob_uncal),
        ("platt", prob_platt),
        ("isotonic", prob_iso),
    ]:
        _, bins_df = expected_calibration_error(y, prob, n_bins=10)
        bins_df.insert(0, "calibration_method", name)
        bin_tables.append(bins_df)

    all_bins = pd.concat(bin_tables, ignore_index=True)
    all_bins.to_csv(OUT_BINS, index=False)

    # Reliability diagram
    plt.figure(figsize=(6.5, 5.5))
    plt.plot([0, 1], [0, 1], linestyle="--", label="Perfect calibration")

    for label, prob in [
        ("Uncalibrated", prob_uncal),
        ("Platt scaling", prob_platt),
        ("Isotonic regression", prob_iso),
    ]:
        frac_pos, mean_pred = calibration_curve(y, prob, n_bins=10, strategy="uniform")
        plt.plot(mean_pred, frac_pos, marker="o", label=label)

    plt.xlabel("Mean predicted probability")
    plt.ylabel("Observed fraction positive")
    plt.title("Reliability diagram for local-context DNABERT-2 probe")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=300)
    plt.close()

    with open(OUT_TXT, "w") as f:
        f.write("DNABERT calibration analysis\n")
        f.write("============================\n\n")
        f.write(f"Input metadata: {DATA_PATH}\n")
        f.write(f"Input embeddings: {EMB_PATH}\n")
        f.write(f"Label column: {label_col}\n")
        f.write(f"Samples: {len(y)}\n")
        f.write(f"Positive labels: {int(y.sum())}\n")
        f.write(f"Negative labels: {int((1-y).sum())}\n\n")
        f.write(metrics.to_string(index=False))
        f.write("\n\nInterpretation note:\n")
        f.write(
            "Calibration metrics are computed from out-of-fold predicted probabilities. "
            "ECE uses 10 equal-width probability bins. Calibration should be interpreted "
            "as probability reliability for weak functional severity labels, not as "
            "experimental safety validation.\n"
        )

    print("\nSaved:")
    print(OUT_TABLE)
    print(OUT_BINS)
    print(OUT_FIG)
    print(OUT_TXT)

    print("\nCalibration metrics:")
    print(metrics.to_string(index=False))

if __name__ == "__main__":
    main()
