# CRISPR-Cas Mediated Re-Sensitization of Multidrug-Resistant Bacteria

**ACS Spring 2026 | BIOT: Division of Biochemical Technology | Paper ID: 4421933**  
**Session:** Microbial Engineering and Fermentation  
**Date:** March 25, 2026  

---

## Overview

This repository contains a computational pipeline for designing, scoring, and ranking **CRISPR-Cas9 guide RNAs (gRNAs)** targeting key antimicrobial resistance (AMR) genes.

The objective is to identify high-specificity, high-efficiency guides capable of **re-sensitizing multidrug-resistant (MDR) bacteria** by disrupting resistance genes at the DNA level.

A total of **6,373 guide RNA candidates** were screened across four critical MDR genes.

---

## Dataset Summary

| Gene | Resistance Type | Total Guides | Excellent Guides | Best Score |
|------|----------------|-------------:|-----------------:|-----------:|
| blaKPC | Carbapenem | 136 | 129 (94.9%) | 90.4 |
| blaNDM1 | Carbapenem | 157 | 153 (97.5%) | 90.4 |
| mcr1 | Colistin (last-resort) | 5,962 | 5,004 (83.9%) | 90.4 |
| mecA | Methicillin / MRSA | 118 | 74 (62.7%) | 90.4 |

---

## Background

Antimicrobial resistance (AMR) causes approximately **1.27 million deaths annually**, with projections reaching **10 million deaths per year by 2050**.

CRISPR-Cas systems can be repurposed to selectively target and cleave resistance genes, offering a **precision alternative to antibiotics**.

This project computationally identifies optimal guide RNAs for four high-priority resistance genes to support future experimental validation.

---

## Pipeline
Gene Sequence (FASTA)
↓
PAM Detection (NGG motifs)
↓
Spacer Extraction (20-nt)
↓
On-Target Scoring
↓
Off-Target Screening (0–3 mismatches)
↓
Conservation Analysis (multi-strain)
↓
Final Composite Scoring
↓
Ranked Guide RNA Candidates

---

## Scoring Model
Final Score = (On-Target × 0.45) + (Specificity × 0.30) + (Conservation × 0.25)

### Specificity Score
Specificity = 100 * exp(-penalty / 50)

Where:
Penalty = (0mm × 50) + (1mm × 20) + (2mm × 8) + (3mm × 3)

---

## Classification

| Classification | Score |
|----------------|------|
| Excellent | ≥ 80 |
| Good | 65 – 79 |
| Moderate | 50 – 64 |
| Poor | < 50 |

---

## Repository Structure
crispr-mdr-resensitization/
│
├── data/
│ ├── targets_multistrain/ # Multi-strain FASTA sequences
│ ├── genomes/ # Background genomes
│ └── plasmids/ # Plasmid sequences
│
├── results_panstrain/
│ ├── all_panstrain_guide_candidates.csv
│ ├── top20_per_gene_panstrain.csv
│ ├── top30_global_panstrain_guides.csv
│ ├── summary_statistics_panstrain.csv
│ ├── Figure1_PanStrainScoreDistribution.png
│ ├── Figure2_TopPanStrainGuides_Global.png
│ ├── Figure3_PanStrainGeneComparison.png
│ ├── Figure4_Specificity_vs_Conservation.png
│ ├── Figure5_GuideCoverageByGene.png
│ └── Heatmap_<gene>_TopGuides.png
│
├── src/
│ └── main_pipeline.py # Full pipeline script
│
├── requirements.txt
└── README.md

---

## Output Columns

| Column | Description |
|--------|-------------|
| gene | Target gene |
| position | Genomic position |
| strand | + / - |
| spacer | 20-nt guide sequence |
| pam | PAM sequence |
| gc_content | GC % |
| on_target_score | Efficiency score |
| offtarget_hits_0mm | Exact matches |
| offtarget_hits_1mm | 1 mismatch hits |
| offtarget_hits_2mm | 2 mismatch hits |
| offtarget_hits_3mm | 3 mismatch hits |
| offtarget_penalty | Weighted penalty |
| specificity_score | Off-target specificity |
| conservation_score | Multi-strain coverage |
| final_score | Composite score |
| classification | Guide quality |

---

## Top Guides (Score: 90.4)

| Gene | Position | Spacer | PAM | GC% |
|------|---------:|--------|-----|----:|
| blaKPC | 756 | CACAATAGGTGCGCGCCCAG | TGG | 65 |
| blaNDM1 | 677 | CGGTGATATTGTCACTGGTG | TGG | 50 |
| mcr1 | 2650 | AACCTGCGCAGGTAACACAG | TGG | 55 |
| mecA | 1965 | AAAGTGGCAGACAAATTGGG | TGG | 45 |

All top guides show **zero off-target hits (0–3 mismatches)** in background genomes.

---

## Known Limitations

- GC filtering not enforced (40–70% recommended range)
- Minor specificity inconsistency in subset of mcr1 guides
- No minimum on-target threshold applied
- Results are computational (no experimental validation yet)

---

## Requirements

- Python >= 3.8
- numpy
- pandas
- matplotlib

Install:
pip install -r requirements.txt

---

## Run Pipeline
python src/main_pipeline.py

Outputs will be saved in:
results_panstrain/

---

## Research Contribution

This work provides a scalable computational framework for **CRISPR-based antimicrobial design**, focusing on:

- Multi-strain conservation
- Genome-wide off-target safety
- High-confidence guide prioritization

---

## Citation

Basu, S.  
CRISPR-Cas mediated re-sensitization of multidrug-resistant bacteria  
ACS Spring 2026 | BIOT Division | Paper ID: 4421933  

---

## Contact

Sanghati Basu  
GitHub: https://github.com/SANGHATI23  

---

## Disclaimer

This repository contains computational predictions only.  
Experimental validation is required before biological or clinical application.# crispr-mdr-resensitization
