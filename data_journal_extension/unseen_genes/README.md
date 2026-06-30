# Unseen AMR Gene FASTA Inputs

This folder will store curated reference coding sequences for prospective unseen AMR genes used in the journal-extension validation.

Planned genes:
- vanA
- tetM
- ermB
- aac6Ib, manuscript label: aac(6')-Ib
- qnrS

Sequence rules:
1. Use coding sequences only, not whole genomes, for the first guide-enumeration pass.
2. Prefer CARD or NCBI RefSeq sequence records.
3. Preserve sequence provenance in the FASTA header.
4. Do not mix unverified downloaded sequences with generated guide outputs.
5. Keep filenames consistent with unseen_gene_reference_plan.csv.

Expected files:
- vanA.fa
- tetM.fa
- ermB.fa
- aac6Ib.fa
- qnrS.fa
