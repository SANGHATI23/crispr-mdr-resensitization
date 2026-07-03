# Reviewer-strengthening experiments

This folder contains two additional computational validation experiments added for the BIBM/FOTR-CRISPR manuscript.

## Experiment 1: CARD/RGI + AMRFinderPlus orthogonal AMR database validation

Folder: `results_reviewer_strengthening/amrfinderplus_validation/`

Final summary:
- AMRFinderPlus detected 17 NCBI AMR-supported hits.
- AMRFinderPlus confirmed 8 of 30 CARD/RGI Perfect/Strict unique feature-level calls.
- Percent confirmation: 26.67%.

Key file:
- `final_card_rgi_amrfinderplus_validation_summary.tsv`

## Experiment 2: CRISPOR top-10 comparator

Folder: `results_reviewer_strengthening/crispor_comparison/`

Final summary:
- blaKPC: overlap 3/10, Jaccard 0.1765
- blaNDM1: overlap 4/10, Jaccard 0.2500
- mcr1: overlap 1/10, Jaccard 0.0526
- mecA: overlap 10/10, Jaccard 1.0000

Key files:
- `crispor_vs_fotr_v2_top10_jaccard_summary.tsv`
- `crispor_vs_fotr_v2_top10_overlap_details.tsv`

Claim boundary:
The CRISPOR comparison is a guide-ranking comparator against CRISPOR output on submitted sequences, not bacterial genome-wide CRISPOR off-target validation.
