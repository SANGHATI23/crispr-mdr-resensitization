# Score-vs-safety divergence analysis

This reviewer-strengthening analysis quantifies whether ranking guides by raw model/composite score gives the same choices as ranking guides by conserved off-target burden.

## Scope

This run compares the current final risk-aware top-5 shortlist for each target gene. It does **not** claim to cover the full 145-candidate pool, because conserved functional off-target burden summaries are currently available for the risk-aware shortlist rather than every candidate guide.

## Inputs

- `results_reviewer_strengthening/comparator_overlap/fotr_guides_with_cfd_rs3_comparator_ranks.csv`
- `results_reviewer_strengthening/variance_ci/conserved_burden_variance_ci_by_guide.csv`
- `results_reviewer_strengthening/variance_ci/guide_accession_burden_matrix.csv`

## Ranking definitions

- Raw-score rank: higher `mlcb_base_score` is better, with `final_gene_rank` used as tie-breaker.
- Safety rank: lower `total_accession_burden_sum` is better, followed by lower `mean_accession_burden`, fewer accessions with nonzero burden, fewer off-target hit rows, higher base score, and then original final rank.

## Outputs

- `results_reviewer_strengthening/score_safety_divergence/score_safety_ranked_detail_within_riskaware_shortlist.csv`
- `results_reviewer_strengthening/score_safety_divergence/score_safety_top3_divergence_by_gene.csv`
- `results_reviewer_strengthening/score_safety_divergence/score_safety_top3_divergence_summary.csv`

## Main reviewer-facing interpretation

This analysis converts the single blaKPC anecdote into a gene-level summary across the evaluated panel, showing how often the safety-aware ranking changes the top guide and how much the raw-score top-3 overlaps with the conserved-burden top-3.
