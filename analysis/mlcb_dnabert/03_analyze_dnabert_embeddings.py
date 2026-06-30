#!/usr/bin/env python3
"""
03_analyze_dnabert_embeddings.py

Purpose
-------
Analyze DNABERT-2 guide embeddings for the MLCB paper.

This script:
1. Loads DNABERT-2 embeddings generated in Step 2.
2. Performs PCA dimensionality reduction.
3. Creates publication-ready PCA figures.
4. Runs lightweight ML probes using frozen DNABERT-2 embeddings.
5. Tests whether DNABERT-2 embeddings can recover FOTR-derived guide priority labels.
6. Writes tables and figures for the MLCB paper.

Inputs
------
results_mlcb/embeddings/dnabert2_guide_embeddings.npy
results_mlcb/embeddings/dnabert2_guide_embedding_metadata.csv

Outputs
-------
results_mlcb/figures/figure_dnabert_pca_by_gene.png
results_mlcb/figures/figure_dnabert_pca_by_score_priority.png
results_mlcb/figures/figure_dnabert_pca_by_base_score.png
results_mlcb/tables/dnabert_embedding_pca_coordinates.csv
results_mlcb/tables/dnabert_probe_performance.csv
results_mlcb/tables/dnabert_nearest_neighbors.csv
results_mlcb/tables/dnabert_analysis_summary.txt
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import make_scorer, f1_score, balanced_accuracy_score, accuracy_score
from sklearn.metrics.pairwise import cosine_similarity


warnings.filterwarnings("ignore")


PROJECT_ROOT = Path(__file__).resolve().parents[2]

EMBEDDING_FILE = PROJECT_ROOT / "results_mlcb" / "embeddings" / "dnabert2_guide_embeddings.npy"
METADATA_FILE = PROJECT_ROOT / "results_mlcb" / "embeddings" / "dnabert2_guide_embedding_metadata.csv"

FIGURE_DIR = PROJECT_ROOT / "results_mlcb" / "figures"
TABLE_DIR = PROJECT_ROOT / "results_mlcb" / "tables"

FIGURE_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PCA_BY_GENE = FIGURE_DIR / "figure_dnabert_pca_by_gene.png"
OUTPUT_PCA_BY_PRIORITY = FIGURE_DIR / "figure_dnabert_pca_by_score_priority.png"
OUTPUT_PCA_BY_SCORE = FIGURE_DIR / "figure_dnabert_pca_by_base_score.png"

OUTPUT_PCA_TABLE = TABLE_DIR / "dnabert_embedding_pca_coordinates.csv"
OUTPUT_PROBE_TABLE = TABLE_DIR / "dnabert_probe_performance.csv"
OUTPUT_NEIGHBOR_TABLE = TABLE_DIR / "dnabert_nearest_neighbors.csv"
OUTPUT_SUMMARY = TABLE_DIR / "dnabert_analysis_summary.txt"


RANDOM_STATE = 42


def load_inputs() -> tuple[np.ndarray, pd.DataFrame]:
    """
    Load embedding matrix and metadata.
    """
    if not EMBEDDING_FILE.exists():
        raise FileNotFoundError(
            f"Missing embeddings: {EMBEDDING_FILE}\n"
            "Run:\n"
            "python analysis/mlcb_dnabert/02_generate_dnabert_embeddings.py"
        )

    if not METADATA_FILE.exists():
        raise FileNotFoundError(
            f"Missing metadata: {METADATA_FILE}\n"
            "Run:\n"
            "python analysis/mlcb_dnabert/02_generate_dnabert_embeddings.py"
        )

    embeddings = np.load(EMBEDDING_FILE)
    metadata = pd.read_csv(METADATA_FILE)

    if embeddings.shape[0] != len(metadata):
        raise ValueError(
            f"Embedding rows ({embeddings.shape[0]}) do not match metadata rows ({len(metadata)})."
        )

    return embeddings, metadata


def run_pca(embeddings: np.ndarray, metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Run PCA on standardized embeddings and return coordinate table.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(embeddings)

    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    coords = pca.fit_transform(X_scaled)

    pca_df = metadata.copy()
    pca_df["PC1"] = coords[:, 0]
    pca_df["PC2"] = coords[:, 1]
    pca_df["PC1_explained_variance_ratio"] = pca.explained_variance_ratio_[0]
    pca_df["PC2_explained_variance_ratio"] = pca.explained_variance_ratio_[1]

    return pca_df


def plot_pca_by_gene(pca_df: pd.DataFrame) -> None:
    """
    PCA scatter plot colored by gene.
    """
    plt.figure(figsize=(8, 6))

    genes = sorted(pca_df["mlcb_gene"].astype(str).unique())

    for gene in genes:
        subset = pca_df[pca_df["mlcb_gene"].astype(str) == gene]
        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=gene,
            alpha=0.8,
            s=45,
        )

    pc1_var = pca_df["PC1_explained_variance_ratio"].iloc[0] * 100
    pc2_var = pca_df["PC2_explained_variance_ratio"].iloc[0] * 100

    plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    plt.title("DNABERT-2 Guide Embedding Space by AMR Gene")
    plt.legend(title="Gene", fontsize=9)
    plt.tight_layout()
    plt.savefig(OUTPUT_PCA_BY_GENE, dpi=300)
    plt.close()


def plot_pca_by_priority(pca_df: pd.DataFrame) -> None:
    """
    PCA scatter plot colored by score-derived high-priority label.
    """
    label_col = "label_high_priority_by_score_top_quartile"

    if label_col not in pca_df.columns:
        print(f"Skipping priority PCA plot. Missing column: {label_col}")
        return

    plt.figure(figsize=(8, 6))

    labels = sorted(pca_df[label_col].dropna().unique())

    for label in labels:
        subset = pca_df[pca_df[label_col] == label]
        label_name = "Top-score quartile" if label == 1 else "Lower-score guides"
        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=label_name,
            alpha=0.8,
            s=45,
        )

    pc1_var = pca_df["PC1_explained_variance_ratio"].iloc[0] * 100
    pc2_var = pca_df["PC2_explained_variance_ratio"].iloc[0] * 100

    plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    plt.title("DNABERT-2 Embedding Space by FOTR-Derived Priority")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(OUTPUT_PCA_BY_PRIORITY, dpi=300)
    plt.close()


def plot_pca_by_score(pca_df: pd.DataFrame) -> None:
    """
    PCA scatter plot colored by continuous base score.
    """
    if "mlcb_base_score" not in pca_df.columns:
        print("Skipping score PCA plot. Missing mlcb_base_score.")
        return

    plt.figure(figsize=(8, 6))

    scatter = plt.scatter(
        pca_df["PC1"],
        pca_df["PC2"],
        c=pca_df["mlcb_base_score"],
        alpha=0.85,
        s=45,
    )

    pc1_var = pca_df["PC1_explained_variance_ratio"].iloc[0] * 100
    pc2_var = pca_df["PC2_explained_variance_ratio"].iloc[0] * 100

    plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    plt.title("DNABERT-2 Guide Embedding Space by FOTR Base Score")
    cbar = plt.colorbar(scatter)
    cbar.set_label("FOTR/base score")
    plt.tight_layout()
    plt.savefig(OUTPUT_PCA_BY_SCORE, dpi=300)
    plt.close()


def run_ml_probes(embeddings: np.ndarray, metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Run lightweight supervised probes.

    Target:
    label_high_priority_by_score_top_quartile

    This is NOT experimental validation.
    It tests whether frozen DNABERT-2 embeddings can recover FOTR-derived ranking labels.
    """
    target_col = "label_high_priority_by_score_top_quartile"

    if target_col not in metadata.columns:
        print(f"Skipping ML probes. Missing target column: {target_col}")
        return pd.DataFrame()

    df = metadata.copy()
    y = df[target_col].astype(int).values

    unique_labels = np.unique(y)
    if len(unique_labels) < 2:
        print("Skipping ML probes. Target has only one class.")
        return pd.DataFrame()

    class_counts = pd.Series(y).value_counts()
    min_class_count = class_counts.min()

    n_splits = min(5, int(min_class_count))
    if n_splits < 2:
        print("Skipping ML probes. Not enough samples in minority class.")
        return pd.DataFrame()

    cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    scoring = {
        "accuracy": make_scorer(accuracy_score),
        "balanced_accuracy": make_scorer(balanced_accuracy_score),
        "macro_f1": make_scorer(f1_score, average="macro"),
    }

    models = {
        "majority_baseline_manual": None,
        "dnabert_logistic_regression_probe": Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RANDOM_STATE,
                    ),
                ),
            ]
        ),
        "dnabert_ridge_classifier_probe": Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                (
                    "clf",
                    RidgeClassifier(
                        class_weight="balanced",
                    ),
                ),
            ]
        ),
    }

    rows = []

    majority_label = pd.Series(y).mode()[0]
    majority_pred = np.repeat(majority_label, len(y))

    rows.append(
        {
            "model": "majority_baseline_manual",
            "target": target_col,
            "cv_folds": n_splits,
            "mean_accuracy": accuracy_score(y, majority_pred),
            "std_accuracy": 0.0,
            "mean_balanced_accuracy": balanced_accuracy_score(y, majority_pred),
            "std_balanced_accuracy": 0.0,
            "mean_macro_f1": f1_score(y, majority_pred, average="macro"),
            "std_macro_f1": 0.0,
            "interpretation": "Majority-class baseline without embedding information.",
        }
    )

    for model_name, model in models.items():
        if model is None:
            continue

        result = cross_validate(
            model,
            embeddings,
            y,
            cv=cv,
            scoring=scoring,
            return_train_score=False,
        )

        rows.append(
            {
                "model": model_name,
                "target": target_col,
                "cv_folds": n_splits,
                "mean_accuracy": np.mean(result["test_accuracy"]),
                "std_accuracy": np.std(result["test_accuracy"]),
                "mean_balanced_accuracy": np.mean(result["test_balanced_accuracy"]),
                "std_balanced_accuracy": np.std(result["test_balanced_accuracy"]),
                "mean_macro_f1": np.mean(result["test_macro_f1"]),
                "std_macro_f1": np.std(result["test_macro_f1"]),
                "interpretation": (
                    "Frozen DNABERT-2 embedding probe testing whether pretrained "
                    "sequence representations recover FOTR-derived guide priority."
                ),
            }
        )

    return pd.DataFrame(rows)


def compute_nearest_neighbors(embeddings: np.ndarray, metadata: pd.DataFrame, top_k: int = 5) -> pd.DataFrame:
    """
    Compute nearest-neighbor guide relationships using cosine similarity.

    This helps analyze whether DNABERT groups guides by gene,
    by score, or by sequence similarity.
    """
    sim = cosine_similarity(embeddings)

    rows = []

    for i in range(sim.shape[0]):
        order = np.argsort(sim[i])[::-1]

        # Remove self
        order = [j for j in order if j != i]

        for rank, j in enumerate(order[:top_k], start=1):
            rows.append(
                {
                    "query_embedding_id": metadata.loc[i, "embedding_id"],
                    "query_gene": metadata.loc[i, "mlcb_gene"],
                    "query_spacer": metadata.loc[i, "mlcb_spacer"],
                    "query_score": metadata.loc[i, "mlcb_base_score"],
                    "neighbor_rank": rank,
                    "neighbor_embedding_id": metadata.loc[j, "embedding_id"],
                    "neighbor_gene": metadata.loc[j, "mlcb_gene"],
                    "neighbor_spacer": metadata.loc[j, "mlcb_spacer"],
                    "neighbor_score": metadata.loc[j, "mlcb_base_score"],
                    "cosine_similarity": sim[i, j],
                    "same_gene": metadata.loc[i, "mlcb_gene"] == metadata.loc[j, "mlcb_gene"],
                }
            )

    return pd.DataFrame(rows)


def write_summary(
    embeddings: np.ndarray,
    metadata: pd.DataFrame,
    pca_df: pd.DataFrame,
    probe_df: pd.DataFrame,
    neighbor_df: pd.DataFrame,
) -> None:
    """
    Write text summary for quick reporting.
    """
    lines = []

    lines.append("DNABERT-2 Embedding Analysis Summary")
    lines.append("====================================")
    lines.append("")
    lines.append(f"Embedding matrix shape: {embeddings.shape}")
    lines.append(f"Metadata rows: {len(metadata)}")
    lines.append(f"Unique genes: {metadata['mlcb_gene'].nunique()}")
    lines.append(f"Unique spacers: {metadata['mlcb_spacer'].nunique()}")
    lines.append("")
    lines.append("Gene distribution:")
    for gene, count in metadata["mlcb_gene"].value_counts().items():
        lines.append(f"- {gene}: {count}")

    lines.append("")
    lines.append("PCA variance explained:")
    lines.append(f"- PC1: {pca_df['PC1_explained_variance_ratio'].iloc[0] * 100:.2f}%")
    lines.append(f"- PC2: {pca_df['PC2_explained_variance_ratio'].iloc[0] * 100:.2f}%")

    if "label_high_priority_by_score_top_quartile" in metadata.columns:
        lines.append("")
        lines.append("FOTR-derived score-priority label distribution:")
        for label, count in metadata["label_high_priority_by_score_top_quartile"].value_counts().items():
            label_name = "top-score quartile" if label == 1 else "lower-score guides"
            lines.append(f"- {label_name}: {count}")

    if not probe_df.empty:
        lines.append("")
        lines.append("Embedding probe performance:")
        for _, row in probe_df.iterrows():
            lines.append(
                f"- {row['model']}: "
                f"balanced accuracy={row['mean_balanced_accuracy']:.3f} "
                f"+/- {row['std_balanced_accuracy']:.3f}, "
                f"macro-F1={row['mean_macro_f1']:.3f} "
                f"+/- {row['std_macro_f1']:.3f}"
            )

    if not neighbor_df.empty:
        same_gene_rate = neighbor_df["same_gene"].mean()
        lines.append("")
        lines.append("Nearest-neighbor structure:")
        lines.append(f"- Fraction of top-k nearest neighbors from same gene: {same_gene_rate:.3f}")

    lines.append("")
    lines.append("Interpretation for MLCB paper:")
    lines.append(
        "These analyses treat DNABERT-2 as a frozen pretrained DNA foundation model. "
        "The supervised probes do not establish experimental CRISPR activity. "
        "They test whether pretrained sequence embeddings contain enough information "
        "to recover FOTR-derived guide-priority labels and to organize AMR guide candidates "
        "by gene or score-related structure."
    )

    with open(OUTPUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main() -> None:
    print("Starting DNABERT-2 embedding analysis...")
    print(f"Project root: {PROJECT_ROOT}")

    embeddings, metadata = load_inputs()

    print(f"Loaded embeddings: {embeddings.shape}")
    print(f"Loaded metadata: {metadata.shape}")

    print("Running PCA...")
    pca_df = run_pca(embeddings, metadata)
    pca_df.to_csv(OUTPUT_PCA_TABLE, index=False)

    print("Creating PCA figures...")
    plot_pca_by_gene(pca_df)
    plot_pca_by_priority(pca_df)
    plot_pca_by_score(pca_df)

    print("Running lightweight ML probes...")
    probe_df = run_ml_probes(embeddings, metadata)
    probe_df.to_csv(OUTPUT_PROBE_TABLE, index=False)

    print("Computing nearest-neighbor structure...")
    neighbor_df = compute_nearest_neighbors(embeddings, metadata, top_k=5)
    neighbor_df.to_csv(OUTPUT_NEIGHBOR_TABLE, index=False)

    print("Writing summary...")
    write_summary(
        embeddings=embeddings,
        metadata=metadata,
        pca_df=pca_df,
        probe_df=probe_df,
        neighbor_df=neighbor_df,
    )

    print("\nWrote figures:")
    print(f"- {OUTPUT_PCA_BY_GENE.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_PCA_BY_PRIORITY.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_PCA_BY_SCORE.relative_to(PROJECT_ROOT)}")

    print("\nWrote tables:")
    print(f"- {OUTPUT_PCA_TABLE.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_PROBE_TABLE.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_NEIGHBOR_TABLE.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    if not probe_df.empty:
        print("\nProbe performance preview:")
        print(probe_df.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()