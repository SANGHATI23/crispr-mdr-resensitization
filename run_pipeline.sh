#!/usr/bin/env bash
set -euo pipefail

mkdir -p outputs models
python scripts/build_feature_matrix.py --config config.yaml
python scripts/train_models.py --config config.yaml
python scripts/score_crispr_targetability.py --config config.yaml
python scripts/make_submission_template.py --config config.yaml

echo "Done. See outputs/ and models/."
