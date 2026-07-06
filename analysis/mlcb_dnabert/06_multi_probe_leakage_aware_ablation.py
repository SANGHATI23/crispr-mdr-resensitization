"""
Multi-probe leakage-aware ablation for local-context DNABERT-2 MLCB analysis.

Purpose:
Evaluate whether the 100-bp local-context DNABERT-2 signal is robust across
multiple lightweight probe families, not only logistic regression.

Inputs:
- results_mlcb/tables/mlcb_guide_dataset_with_local_context_embedding_ids.csv
- results_mlcb/embeddings/dnabert2_local_context_embeddings.npy

Outputs:
- results_mlcb/tables/local_context_dnabert_multi_probe_ablation.csv
- results_mlcb/tables/local_context_dnabert_multi_probe_ablation_summary.txt
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings("ignore")


ROOT = Path(".")
TABLE_DIR = ROOT / "results_mlcb" / "tables"
EMBED_DIR = ROOT / "results_mlcb" / "embeddings"

DATASET_PATH = TABLE_DIR / "mlcb_guide_dataset_with_local_context_embedding_ids.csv"
EMBEDDING_PATH = EMBED_DIR / "dnabert2_local_context_embeddings.npy"

OUT_CSV = TABLE_DIR / "local_context_dnabert_multi_probe_ablation.csv"
OUT_SUMMARY = TABLE_DIR / "local_context_dnabert_multi_probe_ablation_summary.txt"


LABEL_CANDIDATES = [
    "label_high_functional_offtarget_risk",
    "weak_functional_severity_label",
    "label_high_functional_severity",
    "label_high_priority_by_score_top_quartile",
]


LEAKAGE_AWARE_ANNOTATION_FEATURES = [
    "position",
    "gc_content",
    "position_azimuthfile",
    "gc_content_azimuthfile",
    "number_of_reference_matches",
]


def find_label_column(df: pd.DataFrame) -> str:
    for col in LABEL_CANDIDATES:
        if col in df.columns:
            values = pd.Series(df[col]).dropna().unique()
            if len(values) >= 2:
                return col

    binary_cols = []
    for col in df.columns:
        values = set(pd.Series(df[col]).dropna().unique())
        if values.issubset({0, 1, 0.0, 1.0, True, False}) and len(values) == 2:
            binary_cols.append(col)

    if binary_cols:
        print("WARNING: Using first detected binary label column:", binary_cols[0])
        return binary_cols[0]

    raise ValueError(
        "Could not find a usable binary label column. "
        f"Tried: {LABEL_CANDIDATES}"
    )


def get_available_annotation_features(df: pd.DataFrame) -> list:
    features = [c for c in LEAKAGE_AWARE_ANNOTATION_FEATURES if c in df.columns]
    if not features:
        raise ValueError(
            "No leakage-aware annotation features found. Expected one or more of: "
            f"{LEAKAGE_AWARE_ANNOTATION_FEATURES}"
        )
    return features


def make_feature_matrix(
    df: pd.DataFrame,
    embeddings: np.ndarray,
    feature_set: str,
    annotation_features: list,
) -> np.ndarray:
    if feature_set == "annotation_only":
        return df[annotation_features].apply(pd.to_numeric, errors="coerce").values

    if feature_set == "dnabert_only":
        return embeddings

    if feature_set == "hybrid":
        annotation_x = df[annotation_features].apply(pd.to_numeric, errors="coerce").values
        return np.hstack([embeddings, annotation_x])

    raise ValueError(f"Unknown feature_set: {feature_set}")


def build_probes():
    return {
        "logistic_regression": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        solver="liblinear",
                        random_state=42,
                    ),
                ),
            ]
        ),
        "ridge_classifier": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "model",
                    RidgeClassifier(
                        class_weight="balanced",
                        random_state=42,
                    ),
                ),
            ]
        ),
        "knn_k3": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "model",
                    KNeighborsClassifier(
                        n_neighbors=3,
                        weights="distance",
                    ),
                ),
            ]
        ),
        "random_forest_depth3": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=300,
                        max_depth=3,
                        min_samples_leaf=3,
                        class_weight="balanced",
                        random_state=42,
                    ),
                ),
            ]
        ),
    }


def safe_auc(model, x_test, y_test):
    try:
        if hasattr(model, "predict_proba"):
            scores = model.predict_proba(x_test)[:, 1]
        elif hasattr(model, "decision_function"):
            scores = model.decision_function(x_test)
        else:
            return np.nan

        if len(np.unique(y_test)) < 2:
            return np.nan

        return roc_auc_score(y_test, scores)
    except Exception:
        return np.nan


def evaluate_probe(name, model, x, y, cv):
    fold_rows = []

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(x, y), start=1):
        x_train, x_test = x[train_idx], x[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        fitted = clone(model)
        fitted.fit(x_train, y_train)

        preds = fitted.predict(x_test)

        fold_rows.append(
            {
                "probe": name,
                "fold": fold_idx,
                "balanced_accuracy": balanced_accuracy_score(y_test, preds),
                "macro_f1": f1_score(y_test, preds, average="macro"),
                "roc_auc": safe_auc(fitted, x_test, y_test),
            }
        )

    return fold_rows


def summarize(rows):
    df = pd.DataFrame(rows)

    summary = (
        df.groupby(["feature_set", "probe"], as_index=False)
        .agg(
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_std=("balanced_accuracy", "std"),
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            roc_auc_mean=("roc_auc", "mean"),
            roc_auc_std=("roc_auc", "std"),
        )
    )

    return summary


def main():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading dataset:", DATASET_PATH)
    df = pd.read_csv(DATASET_PATH)

    print("Loading embeddings:", EMBEDDING_PATH)
    embeddings = np.load(EMBEDDING_PATH)

    if len(df) != embeddings.shape[0]:
        raise ValueError(
            f"Dataset rows ({len(df)}) do not match embedding rows ({embeddings.shape[0]})."
        )

    label_col = find_label_column(df)
    annotation_features = get_available_annotation_features(df)

    y = pd.to_numeric(df[label_col], errors="coerce").astype(int).values

    print("Dataset shape:", df.shape)
    print("Embedding shape:", embeddings.shape)
    print("Label column:", label_col)
    print("Label counts:", dict(pd.Series(y).value_counts().sort_index()))
    print("Leakage-aware annotation features:", annotation_features)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    probes = build_probes()

    feature_sets = {
        "annotation_only": "Leakage-aware biological annotation features",
        "dnabert_only": "100-bp local-context DNABERT-2 embeddings",
        "hybrid": "100-bp local-context DNABERT-2 + leakage-aware annotation features",
    }

    all_fold_rows = []

    for feature_set, description in feature_sets.items():
        x = make_feature_matrix(df, embeddings, feature_set, annotation_features)

        for probe_name, probe_model in probes.items():
            print(f"Running {feature_set} | {probe_name} | X={x.shape}")
            fold_rows = evaluate_probe(probe_name, probe_model, x, y, cv)

            for row in fold_rows:
                row["feature_set"] = feature_set
                row["feature_description"] = description
                row["label_column"] = label_col
                row["n_guides"] = len(df)
                row["n_features"] = x.shape[1]

            all_fold_rows.extend(fold_rows)

    summary = summarize(all_fold_rows)

    summary["balanced_accuracy"] = summary.apply(
        lambda r: f"{r['balanced_accuracy_mean']:.3f} +/- {r['balanced_accuracy_std']:.3f}",
        axis=1,
    )
    summary["macro_f1"] = summary.apply(
        lambda r: f"{r['macro_f1_mean']:.3f} +/- {r['macro_f1_std']:.3f}",
        axis=1,
    )
    summary["roc_auc"] = summary.apply(
        lambda r: f"{r['roc_auc_mean']:.3f} +/- {r['roc_auc_std']:.3f}"
        if not np.isnan(r["roc_auc_mean"])
        else "NA",
        axis=1,
    )

    summary = summary[
        [
            "feature_set",
            "probe",
            "balanced_accuracy",
            "macro_f1",
            "roc_auc",
            "balanced_accuracy_mean",
            "balanced_accuracy_std",
            "macro_f1_mean",
            "macro_f1_std",
            "roc_auc_mean",
            "roc_auc_std",
        ]
    ]

    summary.to_csv(OUT_CSV, index=False)

    best_rows = summary.sort_values(
        ["feature_set", "balanced_accuracy_mean"],
        ascending=[True, False],
    )

    with open(OUT_SUMMARY, "w") as f:
        f.write("Local-context DNABERT-2 multi-probe leakage-aware ablation\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Dataset: {DATASET_PATH}\n")
        f.write(f"Embeddings: {EMBEDDING_PATH}\n")
        f.write(f"Number of guides: {len(df)}\n")
        f.write(f"Embedding shape: {embeddings.shape[0]} x {embeddings.shape[1]}\n")
        f.write(f"Label column: {label_col}\n")
        f.write(f"Label counts: {dict(pd.Series(y).value_counts().sort_index())}\n")
        f.write(f"Leakage-aware annotation features: {annotation_features}\n\n")
        f.write("Summary table:\n")
        f.write(best_rows.to_string(index=False))
        f.write("\n\nInterpretation guide:\n")
        f.write(
            "If DNABERT-only remains above annotation-only across logistic regression, "
            "ridge, k-NN, and limited-depth random forest, the local-context foundation-model "
            "signal is less likely to be a single-probe artifact.\n"
        )

    print("\nSaved:", OUT_CSV)
    print("Saved:", OUT_SUMMARY)
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
