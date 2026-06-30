#!/usr/bin/env python3
"""
02_generate_dnabert_embeddings.py

Purpose
-------
Generate DNABERT-2 embeddings for the MLCB guide dataset.

Input
-----
results_mlcb/tables/mlcb_guide_dataset.csv

Output
------
results_mlcb/embeddings/dnabert2_guide_embeddings.npy
results_mlcb/embeddings/dnabert2_guide_embedding_metadata.csv
results_mlcb/tables/mlcb_guide_dataset_with_embedding_ids.csv

Notes
-----
This script uses DNABERT-2 as a frozen pretrained DNA foundation model.
It does NOT train DNABERT-2.

Model:
zhihan1996/DNABERT-2-117M

If DNABERT-2 cannot load because of package/environment issues,
the script stops with a clear installation message.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_DATASET = PROJECT_ROOT / "results_mlcb" / "tables" / "mlcb_guide_dataset.csv"

OUTPUT_EMBEDDING_DIR = PROJECT_ROOT / "results_mlcb" / "embeddings"
OUTPUT_EMBEDDING_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_TABLE_DIR = PROJECT_ROOT / "results_mlcb" / "tables"
OUTPUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_EMBEDDINGS = OUTPUT_EMBEDDING_DIR / "dnabert2_guide_embeddings.npy"
OUTPUT_METADATA = OUTPUT_EMBEDDING_DIR / "dnabert2_guide_embedding_metadata.csv"
OUTPUT_DATASET_WITH_IDS = OUTPUT_TABLE_DIR / "mlcb_guide_dataset_with_embedding_ids.csv"

MODEL_NAME = "zhihan1996/DNABERT-2-117M"
MAX_LENGTH = 128
BATCH_SIZE = 8


def load_transformers():
    """
    Import Hugging Face dependencies only when needed.
    """
    try:
        from transformers import AutoTokenizer, AutoModel
    except Exception as exc:
        raise ImportError(
            "\nCould not import transformers.\n\n"
            "Install required packages using:\n"
            "pip install transformers accelerate einops tqdm scikit-learn umap-learn\n\n"
            "Also make sure torch is installed:\n"
            "pip install torch\n"
        ) from exc

    return AutoTokenizer, AutoModel


def get_device() -> torch.device:
    """
    Choose best available device.
    """
    if torch.backends.mps.is_available():
        return torch.device("mps")

    if torch.cuda.is_available():
        return torch.device("cuda")

    return torch.device("cpu")


def mean_pool_last_hidden_state(last_hidden_state: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
    """
    Mean-pool token embeddings using attention mask.
    """
    mask = attention_mask.unsqueeze(-1).expand(last_hidden_state.size()).float()
    summed = torch.sum(last_hidden_state * mask, dim=1)
    counts = torch.clamp(mask.sum(dim=1), min=1e-9)
    return summed / counts


def load_dataset() -> pd.DataFrame:
    """
    Load MLCB guide dataset.
    """
    if not INPUT_DATASET.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_DATASET}\n"
            "Run this first:\n"
            "python analysis/mlcb_dnabert/01_prepare_mlcb_dataset.py"
        )

    df = pd.read_csv(INPUT_DATASET)

    required_cols = ["mlcb_gene", "mlcb_spacer", "dnabert_sequence"]
    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns in input dataset: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    df = df.copy()
    df["embedding_id"] = [f"guide_embedding_{i:04d}" for i in range(len(df))]

    return df


def load_model_and_tokenizer(device: torch.device):
    """
    Load DNABERT-2 tokenizer and model.
    """
    AutoTokenizer, AutoModel = load_transformers()

    print(f"Loading tokenizer: {MODEL_NAME}")
    tokenizer = AutoTokenizer.from_pretrained(
        MODEL_NAME,
        trust_remote_code=True,
    )

    print(f"Loading model: {MODEL_NAME}")
    model = AutoModel.from_pretrained(
        MODEL_NAME,
        trust_remote_code=True,
    )

    model.to(device)
    model.eval()

    return tokenizer, model


@torch.no_grad()
def embed_sequences(
    sequences: list[str],
    tokenizer,
    model,
    device: torch.device,
) -> np.ndarray:
    """
    Generate mean-pooled DNABERT-2 embeddings.
    """
    all_embeddings = []

    for start in tqdm(range(0, len(sequences), BATCH_SIZE), desc="Embedding batches"):
        batch_sequences = sequences[start : start + BATCH_SIZE]

        encoded = tokenizer(
            batch_sequences,
            padding=True,
            truncation=True,
            max_length=MAX_LENGTH,
            return_tensors="pt",
        )

        encoded = {key: value.to(device) for key, value in encoded.items()}

        outputs = model(**encoded)

        if hasattr(outputs, "last_hidden_state"):
            last_hidden = outputs.last_hidden_state
        elif isinstance(outputs, tuple):
            last_hidden = outputs[0]
        else:
            raise RuntimeError("Could not find last_hidden_state in model output.")

        pooled = mean_pool_last_hidden_state(
            last_hidden_state=last_hidden,
            attention_mask=encoded["attention_mask"],
        )

        all_embeddings.append(pooled.detach().cpu().numpy())

    return np.vstack(all_embeddings)


def main() -> None:
    print("Starting DNABERT-2 embedding generation...")
    print(f"Project root: {PROJECT_ROOT}")

    df = load_dataset()
    print(f"Loaded MLCB dataset shape: {df.shape}")

    sequences = df["dnabert_sequence"].astype(str).str.upper().tolist()

    print(f"Number of sequences to embed: {len(sequences)}")
    print(f"Example sequence: {sequences[0]}")

    device = get_device()
    print(f"Using device: {device}")

    tokenizer, model = load_model_and_tokenizer(device)

    embeddings = embed_sequences(
        sequences=sequences,
        tokenizer=tokenizer,
        model=model,
        device=device,
    )

    print(f"Embedding matrix shape: {embeddings.shape}")

    np.save(OUTPUT_EMBEDDINGS, embeddings)

    metadata_cols = [
        "embedding_id",
        "mlcb_gene",
        "mlcb_spacer",
        "mlcb_pam",
        "dnabert_sequence",
        "dnabert_sequence_type",
        "mlcb_base_score",
        "label_high_priority_by_score_top_quartile",
        "weak_functional_severity_label",
        "weak_functional_severity_name",
    ]

    metadata_cols = [col for col in metadata_cols if col in df.columns]

    metadata_df = df[metadata_cols].copy()
    metadata_df.to_csv(OUTPUT_METADATA, index=False)

    df.to_csv(OUTPUT_DATASET_WITH_IDS, index=False)

    print("\nWrote:")
    print(f"- {OUTPUT_EMBEDDINGS.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_METADATA.relative_to(PROJECT_ROOT)}")
    print(f"- {OUTPUT_DATASET_WITH_IDS.relative_to(PROJECT_ROOT)}")

    print("\nDone.")


if __name__ == "__main__":
    main()