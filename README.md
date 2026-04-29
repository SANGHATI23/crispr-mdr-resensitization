# CAMDA CRISPR-AMR Prediction Pipeline

This project converts CABBAGE/CAMDA-style AMR genotype and phenotype CSV files into a reproducible antimicrobial resistance prediction workflow.

## What this pipeline does

1. Builds an assembly-level AMR feature matrix from genotype annotations.
2. Merges phenotype labels for pathogen-drug specific training.
3. Trains one classifier per pathogen-drug combination.
4. Handles class imbalance using class-weighted models.
5. Exports cross-validation metrics and feature importance.
6. Adds a CRISPR targetability layer for mecA, blaKPC, blaNDM, and mcr genes.
7. Produces a CAMDA-style prediction file when test metadata is provided.

## Expected input files

Place the CAMDA/CABBAGE CSV files in `data/`:

- `amr_records.csv`
- `amr_records-2.csv`
- `amr_records-3.csv`
- `amr_records-4.csv`
- `amr_records-5.csv`

The pipeline automatically detects column prefixes:

- `pheno_geno_merged-` for merged phenotype-genotype records
- `genotype-` for genotype-only records
- `phenotype-` for phenotype-only records

## Run

```bash
pip install -r requirements.txt
bash run_pipeline.sh
```

## Main outputs

- `outputs/assembly_feature_matrix.csv`
- `outputs/training_matrix.csv`
- `outputs/model_cv_metrics.csv`
- `outputs/feature_importance.csv`
- `outputs/crispr_targetability_scores.csv`
- `outputs/camda_submission_template.csv`

## Scientific positioning

This is not only a black-box AMR classifier. The distinctive research angle is:

> Predict resistance phenotype and identify biologically actionable CRISPR-resensitization targets among resistance-driving genes.

That means the model produces both:

1. R/S prediction for CAMDA-style evaluation.
2. Gene-level interpretability and CRISPR targetability for blaKPC, blaNDM, mcr, and mecA.

## Important note

If CAMDA provides raw assembly FASTA files rather than precomputed AMR annotations, run AMRFinderPlus, RGI/CARD, or ResFinder first to create the genotype annotation table. This pipeline starts from the annotation table stage.
