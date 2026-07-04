# CARD/RGI confidence-cutoff stress test

This reviewer-strengthening analysis stress-tests the AMR annotation confidence threshold by comparing:

1. Perfect-only
2. Perfect + Strict
3. Perfect + Strict + Loose

## Inputs

- `results_reviewer_strengthening/rgi_card/full_rgi_execution/rgi_card_clean_hit_table.csv`
- `results_reviewer_strengthening/amrfinderplus_validation/amrfinderplus_1161_summary_by_gene.tsv`
- `results_reviewer_strengthening/amrfinderplus_validation/amrfinderplus_1161_summary_by_guide.tsv`

## Outputs

- `results_reviewer_strengthening/confidence_cutoff_stress_test/rgi_cutoff_stress_test_by_gene_long.csv`
- `results_reviewer_strengthening/confidence_cutoff_stress_test/rgi_cutoff_stress_test_by_guide_long.csv`
- `results_reviewer_strengthening/confidence_cutoff_stress_test/rgi_cutoff_stress_test_by_gene_pivot.csv`
- `results_reviewer_strengthening/confidence_cutoff_stress_test/rgi_cutoff_stress_test_by_guide_pivot.csv`
- `results_reviewer_strengthening/confidence_cutoff_stress_test/rgi_cutoff_stress_test_summary.csv`

## Interpretation

This analysis shows whether the AMR annotation signal is dependent on including Loose CARD/RGI hits, or whether key conclusions remain visible under stricter Perfect-only or Perfect+Strict thresholds. AMRFinderPlus summaries are added as an independent database context where available.
