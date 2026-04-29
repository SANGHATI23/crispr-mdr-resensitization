
# STEP 10 Final Official CAMDA AMR Model Rerun

## Purpose
This step reruns clean final models only for the official CAMDA AMR species-antibiotic combinations.

## Official combinations
['staphylococcus_aureus', 'klebsiella_pneumoniae', 'acinetobacter_baumannii', 'escherichia_coli']

## Reproducibility
- Seed: 42
- Train/test split: fixed 80/20 stratified split
- Cross-validation: fixed 5-fold StratifiedKFold
- Feature type: AMR gene presence matrix
- Target: y_resistant

## Output files
- step10_final_train_test_metrics.csv
- step10_final_5fold_cv_metrics.csv
- step10_final_predictions.csv
- final trained full-data models
- step10_final_rerun_config.json
