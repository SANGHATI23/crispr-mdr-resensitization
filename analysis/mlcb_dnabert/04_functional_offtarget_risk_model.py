#!/usr/bin/env python3
"""
04_functional_offtarget_risk_model.py

Purpose
-------
Build the MLCB functional off-target risk layer for FOTR-CRISPR.

This script creates a guide-level functional off-target risk model using:

1. Existing guide-level DNABERT-2 embeddings.
2. Existing FOTR/FOTR-like guide scores.
3. Existing CRISOT specificity/escape information where available.
4. Annotation-derived weak functional severity labels.

Important
---------
This script does NOT claim experimental off-target validation.

It creates an annotation-derived / weakly supervised functional risk layer
for the MLCB paper. The correct paper wording is:

"DNABERT-2 embeddings were used as frozen pretrained DNA representations.
A lightweight probe was trained to recover annotation-derived functional
severity / priority labels, and the predicted functional off-target risk
was used to re-rank guide candidates."

Inputs
------
results_mlcb/tables/mlcb_guide_dataset_with_embedding_ids.csv
results_mlcb/embeddings/dnabert2_guide_embeddings.npy
results_external/crisot/merged/crisot_guide_level_summary.csv  optional

Outputs
-------
results_mlcb/tables/functional_offtarget_model_performance.csv
results_mlcb/tables/functional_offtarget_risk_ranked_guides.csv
results_mlcb/tables/functional_offtarget_disagreement_cases.csv
results_mlcb/tables/functional_offtarget_risk_summary.txt
results_mlcb/figures/figure_functional_risk_ranking_shift.png
results_mlcb/figures/figure_functional_risk_by_gene.png
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate, cross_val_predict
from sklearn.metrics import (
    make_scorer,
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)


warnings.filterwarnings("ignore")


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_GUIDE_DATASET = PROJECT_ROOT / "results_mlcb" / "tables" / "mlcb_guide_dataset_with_embedding_ids.csv"
INPUT_EMBEDDINGS = PROJECT_ROOT / "results_mlcb" / "embeddings" / "dnabert2_guide_embeddings.npy"
INPUT_CRISOT = PROJECT_ROOT / "results_external" / "crisot" / "merged" / "crisot_guide_level_summary.csv"

TABLE_DIR = PROJECT_ROOT / "results_mlcb" / "tables"
FIGURE_DIR = PROJECT_ROOT / "results_mlcb" / "figures"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_MODEL_PERFORMANCE = TABLE_DIR / "functional_offtarget_model_performance.csv"
OUTPUT_RANKED_GUIDES = TABLE_DIR / "functional_offtarget_risk_ranked_guides.csv"
OUTPUT_DISAGREEMENT = TABLE_DIR / "functional_offtarget_disagreement_cases.csv"
OUTPUT_SUMMARY = TABLE_DIR / "functional_offtarget_risk_summary.txt"

OUTPUT_FIG_RANK_SHIFT = FIGURE_DIR / "figure_functional_risk_ranking_shift.png"
OUTPUT_FIG_RISK_BY_GENE = FIGURE_DIR / "figure_functional_risk_by_gene.png"


RANDOM_STATE = 42


def min_max_scale(series: pd.Series, invert: bool = False) -> pd.Series:
    """
    Min-max scale a numeric series to [0, 1].
    If invert=True, high original value becomes low scaled value.
    """
    s = pd.to_numeric(series, errors="coerce")

    if s.notna().sum() == 0:
        scaled = pd.Series(np.zeros(len(s)), index=s.index)
        return scaled

    min_val = s.min()
    max_val = s.max()

    if max_val == min_val:
        scaled = pd.Series(np.zeros(len(s)), index=s.index)
    else:
        scaled = (s - min_val) / (max_val - min_val)

    if invert:
        scaled = 1 - scaled

    return scaled.fillna(0)


def load_inputs() -> tuple[pd.DataFrame, np.ndarray]:
    """
    Load guide dataset and embeddings.
    """
    if not INPUT_GUIDE_DATASET.exists():
        raise FileNotFoundError(
            f"Missing guide dataset: {INPUT_GUIDE_DATASET}\n"
            "Run:\n"
            "python analysis/mlcb_dnabert/02_generate_dnabert_embeddings.py"
        )

    if not INPUT_EMBEDDINGS.exists():
        raise FileNotFoundError(
            f"Missing embeddings: {INPUT_EMBEDDINGS}\n"
            "Run:\n"
            "python analysis/mlcb_dnabert/02_generate_dnabert_embeddings.py"
        )

    guides = pd.read_csv(INPUT_GUIDE_DATASET)
    embeddings = np.load(INPUT_EMBEDDINGS)

    if len(guides) != embeddings.shape[0]:
        raise ValueError(
            f"Guide rows ({len(guides)}) do not match embedding rows ({embeddings.shape[0]})."
        )

    return guides, embeddings


def merge_crisot_if_available(guides: pd.DataFrame) -> pd.DataFrame:
    """
    Merge guide-level CRISOT specificity/escape summary if available.

    CRISOT table is at unique spacer/PAM level. The main 145-row table may not have PAM,
    so we merge by spacer first. If duplicates exist, aggregate by spacer.
    """
    guides = guides.copy()

    if not INPUT_CRISOT.exists():
        guides["crisot_available_for_mlcb"] = False
        guides["crisot_wt_specificity"] = np.nan
        guides["crisot_mutation_count"] = np.nan
        guides["crisot_mean_mutant_score"] = np.nan
        guides["crisot_mean_delta_spec"] = np.nan
        guides["crisot_max_delta_spec"] = np.nan
        return guides

    crisot = pd.read_csv(INPUT_CRISOT)
    crisot.columns = (
        crisot.columns.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(" ", "_", regex=False)
    )

    if "spacer" not in crisot.columns:
        guides["crisot_available_for_mlcb"] = False
        return guides

    keep_cols = [
        "spacer",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_mutant_score",
        "crisot_mean_delta_spec",
        "crisot_max_delta_spec",
    ]

    keep_cols = [col for col in keep_cols if col in crisot.columns]
    crisot = crisot[keep_cols].copy()

    numeric_cols = [col for col in keep_cols if col != "spacer"]
    for col in numeric_cols:
        crisot[col] = pd.to_numeric(crisot[col], errors="coerce")

    crisot_agg = crisot.groupby("spacer", as_index=False)[numeric_cols].mean()

    merged = guides.merge(
        crisot_agg,
        how="left",
        left_on="mlcb_spacer",
        right_on="spacer",
    )

    if "spacer" in merged.columns:
        merged = merged.drop(columns=["spacer"])

    merged["crisot_available_for_mlcb"] = merged["crisot_wt_specificity"].notna()

    return merged


def create_annotation_derived_functional_risk(guides: pd.DataFrame) -> pd.DataFrame:
    """
    Create a weak annotation-derived functional off-target severity target.

    This is currently guide-level, not true off-target-site-level.
    It combines:
    - guide base score
    - CRISOT specificity / mutation risk where available
    - GC extremes
    - gene-level AMR criticality prior

    Label:
    0 = low functional risk
    1 = moderate functional risk
    2 = high functional risk

    For the paper, this must be called weakly supervised / annotation-derived.
    """
    df = guides.copy()

    # Base guide score: high score means better guide, so low risk.
    df["risk_component_low_base_score"] = min_max_scale(df["mlcb_base_score"], invert=True)

    # CRISOT specificity: lower specificity means higher risk.
    if "crisot_wt_specificity" in df.columns:
        df["risk_component_crisot_specificity"] = min_max_scale(
            df["crisot_wt_specificity"],
            invert=True,
        )
    else:
        df["risk_component_crisot_specificity"] = 0.0

    # CRISOT delta specificity: higher change under mutation means higher escape/risk.
    if "crisot_max_delta_spec" in df.columns:
        df["risk_component_crisot_delta"] = min_max_scale(
            df["crisot_max_delta_spec"],
            invert=False,
        )
    else:
        df["risk_component_crisot_delta"] = 0.0

    # GC extremes: very low or very high GC can be less desirable.
    gc = pd.to_numeric(df["mlcb_gc_content"], errors="coerce").fillna(50)
    df["risk_component_gc_extreme"] = (np.abs(gc - 50) / 50).clip(0, 1)

    # Gene-level criticality prior for AMR targets.
    # This is not a biological truth label; it is a reproducible domain prior.
    gene = df["mlcb_gene"].astype(str).str.lower()

    gene_prior_map = {
        "blandm1": 0.90,
        "blakpc": 0.90,
        "meca": 0.85,
        "mcr1": 0.80,
    }

    df["risk_component_gene_amr_prior"] = gene.map(gene_prior_map).fillna(0.75)

    # Combine weak risk components.
    # The weighting is explicit and reproducible.
    df["annotation_derived_functional_risk_score"] = (
        0.35 * df["risk_component_low_base_score"]
        + 0.25 * df["risk_component_crisot_specificity"]
        + 0.15 * df["risk_component_crisot_delta"]
        + 0.10 * df["risk_component_gc_extreme"]
        + 0.15 * df["risk_component_gene_amr_prior"]
    )

    # Convert continuous weak risk into 3 severity classes by tertiles.
    q1 = df["annotation_derived_functional_risk_score"].quantile(1 / 3)
    q2 = df["annotation_derived_functional_risk_score"].quantile(2 / 3)

    df["annotation_derived_functional_severity_label"] = np.select(
        [
            df["annotation_derived_functional_risk_score"] <= q1,
            df["annotation_derived_functional_risk_score"] <= q2,
            df["annotation_derived_functional_risk_score"] > q2,
        ],
        [0, 1, 2],
        default=1,
    )

    df["annotation_derived_functional_severity_name"] = df[
        "annotation_derived_functional_severity_label"
    ].map(
        {
            0: "low",
            1: "moderate",
            2: "high",
        }
    )

    # Binary high-severity target for model probing.
    df["label_high_functional_offtarget_risk"] = np.where(
        df["annotation_derived_functional_severity_label"] == 2,
        1,
        0,
    )

    return df


def create_feature_matrices(df: pd.DataFrame, embeddings: np.ndarray) -> dict[str, np.ndarray]:
    """
    Create feature matrices for model comparison.

    1. sequence_numeric_baseline:
       GC, base score, CRISOT metrics.

    2. dnabert_only:
       frozen 768-dimensional DNABERT-2 embeddings.

    3. dnabert_plus_annotation:
       DNABERT-2 embeddings + annotation/risk numeric features.
    """
    numeric_feature_cols = [
        "mlcb_gc_content",
        "mlcb_base_score",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_mutant_score",
        "crisot_mean_delta_spec",
        "crisot_max_delta_spec",
        "risk_component_low_base_score",
        "risk_component_crisot_specificity",
        "risk_component_crisot_delta",
        "risk_component_gc_extreme",
        "risk_component_gene_amr_prior",
    ]

    available_numeric_cols = [col for col in numeric_feature_cols if col in df.columns]

    numeric_df = df[available_numeric_cols].copy()
    numeric_df = numeric_df.apply(pd.to_numeric, errors="coerce")
    numeric_df = numeric_df.fillna(numeric_df.median(numeric_only=True))
    numeric_df = numeric_df.fillna(0)

    X_numeric = numeric_df.values
    X_dnabert = embeddings
    X_hybrid = np.hstack([X_dnabert, X_numeric])

    return {
        "sequence_and_annotation_numeric_baseline": X_numeric,
        "dnabert_embedding_only": X_dnabert,
        "dnabert_plus_annotation_hybrid": X_hybrid,
    }


def evaluate_functional_risk_models(feature_matrices: dict[str, np.ndarray], y: np.ndarray) -> pd.DataFrame:
    """
    Evaluate lightweight classifiers for high functional off-target risk.
    """
    y = y.astype(int)

    class_counts = pd.Series(y).value_counts()
    min_class_count = int(class_counts.min())

    n_splits = min(5, min_class_count)

    if n_splits < 2 or len(np.unique(y)) < 2:
        return pd.DataFrame(
            [
                {
                    "model": "not_evaluated",
                    "feature_set": "not_evaluated",
                    "target": "label_high_functional_offtarget_risk",
                    "cv_folds": 0,
                    "mean_accuracy": np.nan,
                    "std_accuracy": np.nan,
                    "mean_balanced_accuracy": np.nan,
                    "std_balanced_accuracy": np.nan,
                    "mean_macro_f1": np.nan,
                    "std_macro_f1": np.nan,
                    "mean_roc_auc": np.nan,
                    "std_roc_auc": np.nan,
                    "interpretation": "Not enough target-class variation for cross-validation.",
                }
            ]
        )

    cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    scoring = {
        "accuracy": make_scorer(accuracy_score),
        "balanced_accuracy": make_scorer(balanced_accuracy_score),
        "macro_f1": make_scorer(f1_score, average="macro"),
        "roc_auc": "roc_auc",
    }

    models = {
        "logistic_regression_probe": Pipeline(
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
        "ridge_classifier_probe": Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("clf", RidgeClassifier(class_weight="balanced")),
            ]
        ),
    }

    rows = []

    majority_label = pd.Series(y).mode()[0]
    majority_pred = np.repeat(majority_label, len(y))

    rows.append(
        {
            "model": "majority_baseline_manual",
            "feature_set": "no_features",
            "target": "label_high_functional_offtarget_risk",
            "cv_folds": n_splits,
            "mean_accuracy": accuracy_score(y, majority_pred),
            "std_accuracy": 0.0,
            "mean_balanced_accuracy": balanced_accuracy_score(y, majority_pred),
            "std_balanced_accuracy": 0.0,
            "mean_macro_f1": f1_score(y, majority_pred, average="macro"),
            "std_macro_f1": 0.0,
            "mean_roc_auc": 0.5,
            "std_roc_auc": 0.0,
            "interpretation": "Majority-class baseline.",
        }
    )

    for feature_set_name, X in feature_matrices.items():
        for model_name, model in models.items():
            result = cross_validate(
                model,
                X,
                y,
                cv=cv,
                scoring=scoring,
                return_train_score=False,
                error_score=np.nan,
            )

            rows.append(
                {
                    "model": model_name,
                    "feature_set": feature_set_name,
                    "target": "label_high_functional_offtarget_risk",
                    "cv_folds": n_splits,
                    "mean_accuracy": np.nanmean(result["test_accuracy"]),
                    "std_accuracy": np.nanstd(result["test_accuracy"]),
                    "mean_balanced_accuracy": np.nanmean(result["test_balanced_accuracy"]),
                    "std_balanced_accuracy": np.nanstd(result["test_balanced_accuracy"]),
                    "mean_macro_f1": np.nanmean(result["test_macro_f1"]),
                    "std_macro_f1": np.nanstd(result["test_macro_f1"]),
                    "mean_roc_auc": np.nanmean(result["test_roc_auc"]),
                    "std_roc_auc": np.nanstd(result["test_roc_auc"]),
                    "interpretation": (
                        "Weakly supervised functional off-target risk probe. "
                        "The target is annotation-derived, not experimentally measured."
                    ),
                }
            )

    return pd.DataFrame(rows)


def generate_cross_validated_risk_predictions(
    X: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """
    Generate cross-validated predicted high-functional-risk probabilities.

    Uses logistic regression on the hybrid feature matrix.
    """
    y = y.astype(int)

    class_counts = pd.Series(y).value_counts()
    min_class_count = int(class_counts.min())
    n_splits = min(5, min_class_count)

    if n_splits < 2 or len(np.unique(y)) < 2:
        return np.repeat(y.mean(), len(y))

    cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    model = Pipeline(
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
    )

    proba = cross_val_predict(
        model,
        X,
        y,
        cv=cv,
        method="predict_proba",
    )

    # Probability of high-risk class = class 1
    return proba[:, 1]


def create_final_ranking(df: pd.DataFrame, predicted_risk: np.ndarray) -> pd.DataFrame:
    """
    Create final risk-aware guide ranking.

    Original score: mlcb_base_score
    Predicted functional off-target risk: higher = more risky
    New score: base score penalized by functional risk
    """
    ranked = df.copy()

    ranked["predicted_high_functional_offtarget_risk_probability"] = predicted_risk

    ranked["mlcb_base_score_scaled"] = min_max_scale(ranked["mlcb_base_score"], invert=False)

    ranked["functional_risk_penalized_score"] = (
        0.75 * ranked["mlcb_base_score_scaled"]
        + 0.25 * (1 - ranked["predicted_high_functional_offtarget_risk_probability"])
    )

    ranked["original_rank_by_base_score"] = ranked["mlcb_base_score"].rank(
        ascending=False,
        method="min",
    )

    ranked["risk_aware_rank"] = ranked["functional_risk_penalized_score"].rank(
        ascending=False,
        method="min",
    )

    ranked["rank_shift_original_minus_risk_aware"] = (
        ranked["original_rank_by_base_score"] - ranked["risk_aware_rank"]
    )

    ranked = ranked.sort_values(
        ["risk_aware_rank", "original_rank_by_base_score"],
        ascending=[True, True],
    )

    return ranked


def create_disagreement_cases(ranked: pd.DataFrame, n_cases: int = 15) -> pd.DataFrame:
    """
    Identify important ranking disagreement cases.

    Demoted:
    - looked strong originally but dropped after functional risk penalty.

    Promoted:
    - improved after risk-aware scoring.
    """
    df = ranked.copy()

    demoted = df.sort_values("rank_shift_original_minus_risk_aware", ascending=True).head(n_cases).copy()
    demoted["disagreement_type"] = "demoted_after_functional_risk_penalty"

    promoted = df.sort_values("rank_shift_original_minus_risk_aware", ascending=False).head(n_cases).copy()
    promoted["disagreement_type"] = "promoted_after_low_functional_risk"

    cases = pd.concat([demoted, promoted], axis=0, ignore_index=True)

    keep_cols = [
        "disagreement_type",
        "mlcb_gene",
        "mlcb_spacer",
        "mlcb_pam",
        "mlcb_base_score",
        "mlcb_base_score_scaled",
        "annotation_derived_functional_risk_score",
        "annotation_derived_functional_severity_name",
        "predicted_high_functional_offtarget_risk_probability",
        "functional_risk_penalized_score",
        "original_rank_by_base_score",
        "risk_aware_rank",
        "rank_shift_original_minus_risk_aware",
        "mlcb_gc_content",
        "crisot_available_for_mlcb",
        "crisot_wt_specificity",
        "crisot_mutation_count",
        "crisot_mean_delta_spec",
        "crisot_max_delta_spec",
    ]

    keep_cols = [col for col in keep_cols if col in cases.columns]

    return cases[keep_cols]


def plot_ranking_shift(ranked: pd.DataFrame) -> None:
    """
    Plot original rank vs risk-aware rank.
    """
    plt.figure(figsize=(8, 6))

    plt.scatter(
        ranked["original_rank_by_base_score"],
        ranked["risk_aware_rank"],
        alpha=0.75,
        s=45,
    )

    max_rank = max(
        ranked["original_rank_by_base_score"].max(),
        ranked["risk_aware_rank"].max(),
    )

    plt.plot([1, max_rank], [1, max_rank], linestyle="--", linewidth=1)

    plt.xlabel("Original rank by FOTR/base score")
    plt.ylabel("Risk-aware rank after functional off-target penalty")
    plt.title("Ranking Shift After Functional Off-Target Risk Modeling")
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_RANK_SHIFT, dpi=300)
    plt.close()


def plot_risk_by_gene(ranked: pd.DataFrame) -> None:
    """
    Plot predicted functional risk distribution by gene.
    """
    genes = sorted(ranked["mlcb_gene"].astype(str).unique())
    data = [
        ranked.loc[
            ranked["mlcb_gene"].astype(str) == gene,
            "predicted_high_functional_offtarget_risk_probability",
        ].values
        for gene in genes
    ]

    plt.figure(figsize=(8, 6))
    plt.boxplot(data, labels=genes)
    plt.ylabel("Predicted high functional off-target risk probability")
    plt.xlabel("AMR gene")
    plt.title("Predicted Functional Off-Target Risk by AMR Gene")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_RISK_BY_GENE, dpi=300)
    plt.close()


def write_summary(
    df: pd.DataFrame,
    performance_df: pd.DataFrame,
    ranked: pd.DataFrame,
    disagreement: pd.DataFrame,
) -> None:
    """
    Write a paper-ready summary.
    """
    lines = []

    lines.append("Functional Off-Target Risk Modeling Summary")
    lines.append("==========================================")
    lines.append("")
    lines.append(f"Total guides analyzed: {len(df)}")
    lines.append(f"Unique genes: {df['mlcb_gene'].nunique()}")
    lines.append(f"Unique spacers: {df['mlcb_spacer'].nunique()}")
    lines.append(f"CRISOT available rows: {int(df['crisot_available_for_mlcb'].sum())}")
    lines.append("")

    lines.append("Annotation-derived functional severity distribution:")
    for label, count in df["annotation_derived_functional_severity_name"].value_counts().items():
        lines.append(f"- {label}: {count}")

    lines.append("")
    lines.append("High functional off-target risk label distribution:")
    for label, count in df["label_high_functional_offtarget_risk"].value_counts().items():
        label_name = "high risk" if label == 1 else "not high risk"
        lines.append(f"- {label_name}: {count}")

    lines.append("")
    lines.append("Model performance:")
    for _, row in performance_df.iterrows():
        lines.append(
            f"- {row['feature_set']} / {row['model']}: "
            f"balanced accuracy={row['mean_balanced_accuracy']:.3f} "
            f"+/- {row['std_balanced_accuracy']:.3f}, "
            f"macro-F1={row['mean_macro_f1']:.3f} "
            f"+/- {row['std_macro_f1']:.3f}, "
            f"ROC-AUC={row['mean_roc_auc']:.3f} "
            f"+/- {row['std_roc_auc']:.3f}"
        )

    top = ranked.sort_values("risk_aware_rank").head(10)

    lines.append("")
    lines.append("Top 10 risk-aware guide candidates:")
    for _, row in top.iterrows():
        lines.append(
            f"- {row['mlcb_gene']} | {row['mlcb_spacer']} | "
            f"base_score={row['mlcb_base_score']:.3f} | "
            f"predicted_risk={row['predicted_high_functional_offtarget_risk_probability']:.3f} | "
            f"risk_aware_score={row['functional_risk_penalized_score']:.3f}"
        )

    lines.append("")
    lines.append("Interpretation for MLCB paper:")
    lines.append(
        "This analysis adds a functional off-target consequence layer to FOTR-CRISPR. "
        "The labels are annotation-derived weak labels, not experimental measurements. "
        "The ML contribution is the use of frozen DNABERT-2 representations, alone and in "
        "combination with biological annotation features, to estimate functional risk and "
        "re-rank AMR guide candidates under safety-aware constraints."
    )

    lines.append("")
    lines.append("Safe wording:")
    lines.append(
        "Use: 'annotation-derived functional severity', 'weakly supervised risk model', "
        "'foundation-model-augmented risk-aware ranking'. "
        "Avoid: 'experimentally validated off-target prediction' or 'therapeutic safety proven'."
    )

    with open(OUTPUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main() -> None:
    print("Starting functional off-target risk modeling...")
    print(f"Project root: {PROJECT_ROOT}")

    guides, embeddings = load_inputs()
    print(f"Loaded guide dataset: {guides.shape}")
    print(f"Loaded embeddings: {embeddings.shape}")

    guides = merge_crisot_if_available(guides)
    print("Merged CRISOT if available.")
    print(f"CRISOT available rows: {int(guides['crisot_available_for_mlcb'].sum())}")

    guides = create_annotation_derived_functional_risk(guides)
    print("Created annotation-derived functional severity labels.")

    y = guides["label_high_functional_offtarget_risk"].astype(int).values
    feature_matrices = create_feature_matrices(guides, embeddings)

    print("Evaluating functional risk models...")
    performance_df = evaluate_functional_risk_models(feature_matrices, y)
    performance_df.to_csv(OUTPUT_MODEL_PERFORMANCE, index=False)

    print("Generating cross-validated hybrid risk predictions...")
    hybrid_X = feature_matrices["dnabert_plus_annotation_hybrid"]
    predicted_risk = generate_cross_validated_risk_predictions(hybrid_X, y)

    ranked = create_final_ranking(guides, predicted_risk)
    ranked.to_csv(OUTPUT_RANKED_GUIDES, index=False)

    disagreement = create_disagreement_cases(ranked, n_cases=15)
    disagreement.to_csv(OUTPUT_DISAGREEMENT, index=False)

    print("Creating figures...")
    plot_ranking_shift(ranked)
    plot_risk_by_gene(ranked)

    print("Writing summary...")
    write_summary(
        df=guides,
        performance_df=performance_df,
        ranked=ranked,
        disagreement=disagreement,
    )

    print("\nWrote tables:")
    print(f"- {OUTPUT_MODEL_PERFORMANCE.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_RANKED_GUIDES.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_DISAGREEMENT.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_SUMMARY.relative_to(PROJECT_ROOT)}")

    print("\nWrote figures:")
    print(f"- {OUTPUT_FIG_RANK_SHIFT.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_FIG_RISK_BY_GENE.relative_to(PROJECT_ROOT)}")

    print("\nPerformance preview:")
    print(performance_df.to_string(index=False))

    print("\nTop 10 risk-aware guides:")
    preview_cols = [
        "mlcb_gene",
        "mlcb_spacer",
        "mlcb_base_score",
        "annotation_derived_functional_risk_score",
        "predicted_high_functional_offtarget_risk_probability",
        "functional_risk_penalized_score",
        "risk_aware_rank",
    ]
    preview_cols = [col for col in preview_cols if col in ranked.columns]
    print(ranked[preview_cols].head(10).to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()