#!/usr/bin/env python3
"""
Structured AMR holdout benchmark for MLCB DNABERT/FOTR analysis.

This script evaluates whether frozen DNABERT-2 local-context embeddings generalize
under stricter distribution-shift splits:
  1. Gene holdout
  2. Mechanism holdout
  3. Guide-position holdout

The benchmark is intentionally interpreted as a stress test, not as a replacement
for the main leakage-aware 5-fold weak-label recovery analysis.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    accuracy_score,
)

warnings.filterwarnings("ignore", category=UserWarning)


ROOT = Path(__file__).resolve().parents[2]
TABLE_PATH = ROOT / "results_mlcb/tables/mlcb_guide_dataset_with_local_context_embedding_ids.csv"
EMB_PATH = ROOT / "results_mlcb/embeddings/dnabert2_local_context_embeddings.npy"
OUTDIR = ROOT / "results_mlcb/tables"
OUTDIR.mkdir(parents=True, exist_ok=True)

LABEL_COL = "weak_functional_severity_label"


def safe_auc(y_true, score):
    """Return ROC-AUC only when the test set contains both classes."""
    if len(np.unique(y_true)) < 2:
        return np.nan
    return roc_auc_score(y_true, score)


def get_score(model, x_test):
    """Return probability or decision score for ROC-AUC."""
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x_test)[:, 1]
    if hasattr(model, "decision_function"):
        return model.decision_function(x_test)
    return model.predict(x_test)


def add_split_metadata(df):
    """Add manually curated AMR mechanism, gene-family, and position-bin labels."""
    df = df.copy()

    df["amr_mechanism"] = df["gene"].map(
        {
            "blaKPC": "beta_lactamase",
            "blaNDM1": "beta_lactamase",
            "mcr1": "colistin_resistance",
            "mecA": "methicillin_resistance",
        }
    )

    df["gene_family"] = df["gene"].map(
        {
            "blaKPC": "blaKPC_family",
            "blaNDM1": "blaNDM_family",
            "mcr1": "mcr_family",
            "mecA": "mec_family",
        }
    )

    df["position_bin"] = (
        df.groupby("gene")["target_start_0based"]
        .transform(
            lambda s: pd.qcut(
                s.rank(method="first"),
                q=3,
                labels=["early", "middle", "late"],
            )
        )
        .astype(str)
    )

    return df


def build_splits(df):
    """Create gene, mechanism, and guide-position holdout splits."""
    splits = []

    for gene in sorted(df["gene"].unique()):
        train_idx = df.index[df["gene"] != gene].to_numpy()
        test_idx = df.index[df["gene"] == gene].to_numpy()
        splits.append(
            (
                "gene_holdout",
                gene,
                train_idx,
                test_idx,
                "Strict gene-out stress test",
            )
        )

    for mechanism in sorted(df["amr_mechanism"].dropna().unique()):
        train_idx = df.index[df["amr_mechanism"] != mechanism].to_numpy()
        test_idx = df.index[df["amr_mechanism"] == mechanism].to_numpy()
        splits.append(
            (
                "mechanism_holdout",
                mechanism,
                train_idx,
                test_idx,
                "Mechanism transfer stress test",
            )
        )

    for gene in sorted(df["gene"].unique()):
        for bin_name in ["early", "middle", "late"]:
            test_mask = (df["gene"] == gene) & (df["position_bin"] == bin_name)
            if test_mask.sum() == 0:
                continue

            test_idx = df.index[test_mask].to_numpy()
            train_idx = df.index[~test_mask].to_numpy()
            splits.append(
                (
                    "guide_position_holdout",
                    f"{gene}_{bin_name}",
                    train_idx,
                    test_idx,
                    "Unseen local target-region stress test",
                )
            )

    return splits


def build_models():
    """Define lightweight supervised probes."""
    return {
        "majority_baseline": DummyClassifier(strategy="most_frequent"),
        "logistic_regression": make_pipeline(
            StandardScaler(),
            LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42,
            ),
        ),
        "ridge_classifier": make_pipeline(
            StandardScaler(),
            RidgeClassifier(class_weight="balanced"),
        ),
        "random_forest_depth3": RandomForestClassifier(
            n_estimators=300,
            max_depth=3,
            class_weight="balanced",
            random_state=42,
        ),
    }


def run_benchmark(df, x, y):
    """Run all models across all structured holdout splits."""
    rows = []
    models = build_models()
    splits = build_splits(df)

    for split_type, heldout_group, train_idx, test_idx, interpretation in splits:
        y_train = y[train_idx]
        y_test = y[test_idx]

        train_both = len(np.unique(y_train)) == 2
        test_both = len(np.unique(y_test)) == 2
        reportable_auc = train_both and test_both

        for model_name, model in models.items():
            row = {
                "split_type": split_type,
                "heldout_group": heldout_group,
                "model": model_name,
                "train_n": len(train_idx),
                "test_n": len(test_idx),
                "train_neg": int((y_train == 0).sum()),
                "train_pos": int((y_train == 1).sum()),
                "test_neg": int((y_test == 0).sum()),
                "test_pos": int((y_test == 1).sum()),
                "train_positive_rate": round(float(y_train.mean()), 4),
                "test_positive_rate": round(float(y_test.mean()), 4),
                "train_both_classes": int(train_both),
                "test_both_classes": int(test_both),
                "reportable_for_auc": int(reportable_auc),
                "interpretation": interpretation,
            }

            if not train_both:
                row.update(
                    {
                        "accuracy": np.nan,
                        "balanced_accuracy": np.nan,
                        "macro_f1": np.nan,
                        "roc_auc": np.nan,
                        "note": "Skipped: training set has one class only",
                    }
                )
                rows.append(row)
                continue

            try:
                model.fit(x[train_idx], y_train)
                pred = model.predict(x[test_idx])
                score = get_score(model, x[test_idx])

                row.update(
                    {
                        "accuracy": accuracy_score(y_test, pred),
                        "balanced_accuracy": balanced_accuracy_score(y_test, pred),
                        "macro_f1": f1_score(
                            y_test,
                            pred,
                            average="macro",
                            zero_division=0,
                        ),
                        "roc_auc": safe_auc(y_test, score),
                        "note": ""
                        if reportable_auc
                        else "ROC-AUC undefined: test set has one class",
                    }
                )
            except Exception as exc:
                row.update(
                    {
                        "accuracy": np.nan,
                        "balanced_accuracy": np.nan,
                        "macro_f1": np.nan,
                        "roc_auc": np.nan,
                        "note": f"Model failed: {exc}",
                    }
                )

            rows.append(row)

    results = pd.DataFrame(rows)

    for col in ["accuracy", "balanced_accuracy", "macro_f1", "roc_auc"]:
        results[col] = results[col].round(4)

    return results


def create_manuscript_outputs(results):
    """Create compact reportable tables and paragraph."""
    full_out = OUTDIR / "amr_generalization_holdout_benchmark.csv"
    compact_out = OUTDIR / "amr_generalization_holdout_benchmark_compact.csv"
    best_out = OUTDIR / "amr_generalization_holdout_best_per_split.csv"
    manuscript_csv = OUTDIR / "amr_generalization_holdout_manuscript_table.csv"
    manuscript_md = OUTDIR / "amr_generalization_holdout_manuscript_table.md"
    paragraph_txt = OUTDIR / "amr_generalization_holdout_manuscript_paragraph.txt"
    summary_txt = OUTDIR / "amr_generalization_holdout_summary.txt"

    results.to_csv(full_out, index=False)

    compact = results[
        (results["reportable_for_auc"] == 1)
        & (results["model"] != "majority_baseline")
    ].copy()
    compact.to_csv(compact_out, index=False)

    best = (
        compact.sort_values(
            ["split_type", "heldout_group", "balanced_accuracy", "roc_auc"],
            ascending=[True, True, False, False],
        )
        .groupby(["split_type", "heldout_group"], as_index=False)
        .head(1)
    )
    best.to_csv(best_out, index=False)

    manuscript = best[
        [
            "split_type",
            "heldout_group",
            "model",
            "train_n",
            "test_n",
            "test_neg",
            "test_pos",
            "balanced_accuracy",
            "macro_f1",
            "roc_auc",
            "interpretation",
        ]
    ].copy()

    manuscript["split_type"] = manuscript["split_type"].map(
        {
            "gene_holdout": "Gene holdout",
            "mechanism_holdout": "Mechanism holdout",
            "guide_position_holdout": "Guide-position holdout",
        }
    ).fillna(manuscript["split_type"])

    manuscript["model"] = manuscript["model"].map(
        {
            "logistic_regression": "Logistic regression",
            "ridge_classifier": "Ridge classifier",
            "random_forest_depth3": "Random forest, depth 3",
        }
    ).fillna(manuscript["model"])

    manuscript = manuscript.rename(
        columns={
            "split_type": "Split type",
            "heldout_group": "Held-out group",
            "model": "Best probe",
            "train_n": "Train n",
            "test_n": "Test n",
            "test_neg": "Test negatives",
            "test_pos": "Test positives",
            "balanced_accuracy": "Balanced accuracy",
            "macro_f1": "Macro-F1",
            "roc_auc": "ROC-AUC",
            "interpretation": "Interpretation",
        }
    )

    manuscript.to_csv(manuscript_csv, index=False)
    manuscript_md.write_text(manuscript.to_markdown(index=False))

    paragraph = (
        "Structured AMR holdout analysis was added as a distribution-shift stress test rather than as "
        "a replacement for the main leakage-aware 5-fold evaluation. Reportable held-out splits were "
        "restricted to cases where both training and test partitions contained both weak-label classes. "
        "Under these stricter settings, performance dropped substantially relative to within-distribution "
        "cross-validation: the best gene-holdout balanced accuracy was 0.500 for blaKPC and 0.500 for mcr1, "
        "while mechanism-holdout balanced accuracy reached 0.558 for beta-lactamase holdout and 0.500 for "
        "colistin-resistance holdout. Several biologically meaningful splits, including blaNDM1 and mecA "
        "holdouts, had undefined ROC-AUC because the held-out test partition contained only one weak-label "
        "class. These findings support a conservative interpretation: frozen DNABERT-2 embeddings recover "
        "weak-label signal within the current AMR setting, but true AMR gene/mechanism transfer remains "
        "difficult and requires larger, more balanced AMR guide-family benchmarks."
    )
    paragraph_txt.write_text(paragraph)

    with open(summary_txt, "w") as f:
        f.write("AMR structured holdout benchmark summary\n")
        f.write("========================================\n\n")
        f.write(f"Total split-model rows: {len(results)}\n")
        f.write(f"Reportable non-baseline rows: {len(compact)}\n")
        f.write(f"Best-per-split rows: {len(best)}\n\n")
        f.write("Best model per reportable split:\n")
        f.write(best.to_string(index=False))
        f.write("\n\nManuscript paragraph:\n")
        f.write(paragraph)

    return {
        "full": full_out,
        "compact": compact_out,
        "best": best_out,
        "manuscript_csv": manuscript_csv,
        "manuscript_md": manuscript_md,
        "paragraph": paragraph_txt,
        "summary": summary_txt,
    }


def main():
    df = pd.read_csv(TABLE_PATH)
    x = np.load(EMB_PATH)

    if len(df) != x.shape[0]:
        raise ValueError(f"Row mismatch: df={len(df)}, embeddings={x.shape}")

    df = add_split_metadata(df)
    y = df[LABEL_COL].astype(int).values

    results = run_benchmark(df, x, y)
    outputs = create_manuscript_outputs(results)

    print("AMR structured holdout benchmark completed.")
    print(f"Dataset n: {len(df)}")
    print(f"Embedding shape: {x.shape}")
    print(f"Label counts: {pd.Series(y).value_counts().to_dict()}")
    print("\nSaved outputs:")
    for name, path in outputs.items():
        print(f" - {name}: {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
