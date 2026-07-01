#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score, accuracy_score

DATA_PATH = "results_mlcb/tables/mlcb_guide_dataset_with_local_context_embedding_ids.csv"
EMB_PATH = "results_mlcb/embeddings/dnabert2_local_context_embeddings.npy"

OUT_DETAIL = "results_mlcb/tables/cross_gene_heldout_dnabert_results.csv"
OUT_SUMMARY = "results_mlcb/tables/cross_gene_heldout_dnabert_summary.csv"
OUT_TXT = "results_mlcb/tables/cross_gene_heldout_dnabert_summary.txt"

LABEL_CANDIDATES = [
    "weak_functional_severity_label",
    "label_high_functional_offtarget_risk",
    "label_high_priority_by_score_top_quartile",
]

GENE_CANDIDATES = [
    "mlcb_gene",
    "gene",
    "target_gene",
]

def find_col(df, candidates, coltype):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Could not find {coltype} column. Tried: {candidates}. Available columns: {list(df.columns)}")

def safe_auc(y_true, y_score):
    try:
        if len(np.unique(y_true)) < 2:
            return np.nan
        return roc_auc_score(y_true, y_score)
    except Exception:
        return np.nan

def main():
    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(f"Missing metadata table: {DATA_PATH}")
    if not os.path.exists(EMB_PATH):
        raise FileNotFoundError(f"Missing embedding matrix: {EMB_PATH}")

    df = pd.read_csv(DATA_PATH)
    X = np.load(EMB_PATH)

    if len(df) != X.shape[0]:
        raise ValueError(f"Row mismatch: metadata has {len(df)} rows but embeddings have {X.shape[0]} rows")

    gene_col = find_col(df, GENE_CANDIDATES, "gene")
    label_col = find_col(df, LABEL_CANDIDATES, "label")

    y = df[label_col].astype(int).values
    genes = df[gene_col].astype(str).values

    print("Loaded metadata:", df.shape)
    print("Loaded embeddings:", X.shape)
    print("Gene column:", gene_col)
    print("Label column:", label_col)
    print("Gene counts:")
    print(df[gene_col].value_counts())
    print("Label counts:")
    print(df[label_col].value_counts())

    models = {
        "Majority baseline": None,
        "Logistic regression": make_pipeline(
            StandardScaler(),
            LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42,
            )
        ),
        "Ridge classifier": make_pipeline(
            StandardScaler(),
            RidgeClassifier(class_weight="balanced")
        ),
        "k-NN, k=3": make_pipeline(
            StandardScaler(),
            KNeighborsClassifier(n_neighbors=3)
        ),
        "Random forest, depth=3": RandomForestClassifier(
            n_estimators=300,
            max_depth=3,
            class_weight="balanced",
            random_state=42,
        ),
    }

    rows = []

    for heldout_gene in sorted(np.unique(genes)):
        train_idx = genes != heldout_gene
        test_idx = genes == heldout_gene

        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        print(f"\nHeld-out gene: {heldout_gene}")
        print("Train size:", len(y_train), "Test size:", len(y_test))
        print("Train labels:", dict(zip(*np.unique(y_train, return_counts=True))))
        print("Test labels:", dict(zip(*np.unique(y_test, return_counts=True))))

        for model_name, model in models.items():
            if model is None:
                majority_class = pd.Series(y_train).mode().iloc[0]
                y_pred = np.full_like(y_test, majority_class)
                y_score = np.full_like(y_test, float(majority_class), dtype=float)
            else:
                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)

                if hasattr(model, "predict_proba"):
                    y_score = model.predict_proba(X_test)[:, 1]
                elif hasattr(model, "decision_function"):
                    y_score = model.decision_function(X_test)
                else:
                    y_score = y_pred.astype(float)

            rows.append({
                "heldout_gene": heldout_gene,
                "model": model_name,
                "n_train": len(y_train),
                "n_test": len(y_test),
                "test_positive_rate": float(np.mean(y_test)),
                "accuracy": accuracy_score(y_test, y_pred),
                "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
                "macro_f1": f1_score(y_test, y_pred, average="macro", zero_division=0),
                "roc_auc": safe_auc(y_test, y_score),
            })

    results = pd.DataFrame(rows)
    results.to_csv(OUT_DETAIL, index=False)

    summary = (
        results
        .groupby("model", as_index=False)
        .agg(
            mean_balanced_accuracy=("balanced_accuracy", "mean"),
            std_balanced_accuracy=("balanced_accuracy", "std"),
            mean_macro_f1=("macro_f1", "mean"),
            std_macro_f1=("macro_f1", "std"),
            mean_roc_auc=("roc_auc", "mean"),
            std_roc_auc=("roc_auc", "std"),
        )
        .sort_values("mean_balanced_accuracy", ascending=False)
    )

    summary.to_csv(OUT_SUMMARY, index=False)

    with open(OUT_TXT, "w") as f:
        f.write("Cross-gene held-out DNABERT evaluation\n")
        f.write("======================================\n\n")
        f.write(f"Input metadata: {DATA_PATH}\n")
        f.write(f"Input embeddings: {EMB_PATH}\n")
        f.write(f"Gene column: {gene_col}\n")
        f.write(f"Label column: {label_col}\n")
        f.write(f"Embedding shape: {X.shape}\n\n")
        f.write("Per-model summary across held-out genes:\n")
        f.write(summary.to_string(index=False))
        f.write("\n\nInterpretation note:\n")
        f.write(
            "This leave-one-gene-out analysis tests whether DNABERT-2 local-context "
            "embeddings generalize across AMR gene boundaries rather than only learning "
            "within-gene patterns. ROC-AUC is reported as NaN for held-out genes where "
            "the test split contains only one class.\n"
        )

    print("\nSaved:")
    print(OUT_DETAIL)
    print(OUT_SUMMARY)
    print(OUT_TXT)
    print("\nSummary:")
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
