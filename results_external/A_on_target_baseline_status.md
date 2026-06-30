# A. On-target-only baseline completion status

## Completed
1. Internal on-target scoring is complete.
2. Rule Set 3 external-style on-target scoring is complete.
3. RS3 score, internal on-target score, and rank comparison are available.
4. Master cross-model table includes RS3, CFD, conservation, final score, and ranks.

## Limitation
The RS3/Azimuth input used standardized 30mer context rather than perfect true genomic flanking context. Therefore this is a limited-context external benchmark, not a perfect true-context Doench/Azimuth benchmark.

## Not completed
1. Doench 2016 / Azimuth legacy model was not completed because of package/version incompatibility.
2. Perfect true-context RS3 was not completed unless guide flanking sequence is reconstructed directly from the original reference FASTA.

## Final interpretation
A is substantially complete for publication-style benchmarking, but should be reported transparently as:
"Internal on-target baseline plus limited-context Rule Set 3 external comparator."

## Updated completion
A. On-target-only baseline: 90%
