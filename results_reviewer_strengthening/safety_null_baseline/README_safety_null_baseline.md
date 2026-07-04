# Safety null baseline

This analysis tests whether the risk-aware safety layer performs better than a matched random baseline.

## Design

For each target gene, five non-selected candidate guides were randomly sampled from the same candidate guide pool used by the main pipeline. The selected risk-aware top-5 guides were excluded before sampling.

The null guides were screened with Cas-OFFinder against the same 36-accession expanded genome panel and the same PAM/mismatch configuration used for the selected-guide expanded-panel run.

## Inputs

- `results/cross_model_reference_ranked_guides.csv`
- `results_final_selection/cas_offinder_input_guides_top5_per_gene.tsv`
- `results_reviewer_strengthening/safety_null_baseline/null_casoffinder_run_manifest.csv`
- `results_reviewer_strengthening/safety_null_baseline/null_casoffinder_query_mapping.tsv`
- `results_cas_offinder/expanded_panel/offtarget_conservation_requirement4/requirement4_guide_level_offtarget_conservation_ALL20.csv`

## Outputs

- `results_reviewer_strengthening/safety_null_baseline/null_baseline_parsed_offtarget_hits.csv`
- `results_reviewer_strengthening/safety_null_baseline/null_baseline_mapped_offtarget_hits.csv`
- `results_reviewer_strengthening/safety_null_baseline/null_baseline_guide_level_summary.csv`
- `results_reviewer_strengthening/safety_null_baseline/selected_vs_null_safety_burden_by_gene.csv`
- `results_reviewer_strengthening/safety_null_baseline/selected_vs_null_safety_burden_overall_summary.csv`

## Important scope note

This first null baseline compares genomewide Cas-OFFinder hit burden and a transparent mismatch-weighted burden proxy. It does not yet rerun the full database-refined functional annotation layer for the null guides. Therefore, it is appropriate as a direct answer to whether clean/low-hit guides could have appeared by chance, while database-refined null annotation can be added as a further extension.
