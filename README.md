# CRISPR-Cas Mediated Re-Sensitization of Multidrug-Resistant Bacteria

> **ACS Spring 2026** | BIOT: Division of Biochemical Technology | Paper ID: 4421933
> **Session:** Microbial Engineering and Fermentation | March 25, 2026
> **Author:** Sanghati Basu | University of Illinois Springfield

---

## Overview

This repository presents a computational pipeline for designing, scoring, ranking, and validating CRISPR-Cas9 guide RNA (gRNA) candidates targeting major antimicrobial resistance (AMR) genes in multidrug-resistant (MDR) bacteria.

The project focuses on four clinically critical resistance genes:

| Gene | Resistance Type | Clinical Significance |
|------|----------------|----------------------|
| **blaKPC** | Carbapenem | Last-line antibiotic resistance in Gram-negative pathogens |
| **blaNDM1** | Carbapenem | New Delhi metallo-beta-lactamase, globally disseminated |
| **mcr1** | Colistin | Plasmid-mediated last-resort antibiotic resistance |
| **mecA** | Methicillin / MRSA | Major cause of healthcare-associated infections |

The pipeline proceeds in two phases:
1. **Initial computational screen** — broad candidate generation across single reference sequences
2. **Pan-strain validation** — conservation-aware scoring across multiple clinical sequence variants

---

## Results Summary

### Phase 1 — Initial Screen

A total of **6,373 guide RNA candidates** were identified and scored across all four target genes using single-reference genome sequences.

| Gene | Resistance | Total | Excellent | Good | Best Score | Mean Score |
|------|-----------|-------|-----------|------|------------|------------|
| blaKPC | Carbapenem | 136 | 129 (94.9%) | 7 | 90.4 | 83.6 |
| blaNDM1 | Carbapenem | 157 | 153 (97.5%) | 3 | 90.4 | 83.9 |
| mcr1 | Colistin | 5,962 | 5,004 (83.9%) | 750 | 90.4 | 79.6 |
| mecA | Methicillin/MRSA | 118 | 74 (62.7%) | 35 | 90.4 | 72.8 |

### Phase 2 — Pan-Strain Validation

Following strict GC content filtering and multi-strain conservation analysis, **145 elite guides** were retained — each validated across 3 clinical sequence variants per gene.

| Gene | Resistance | Total | Excellent | Good | Best Score | Mean Score | Mean Conservation |
|------|-----------|-------|-----------|------|------------|------------|-------------------|
| blaKPC | Carbapenem | 34 | 34 (100%) | 0 | 86.9 | 84.9 | 80.0 |
| blaNDM1 | Carbapenem | 46 | 42 (91.3%) | 4 | 87.8 | 84.8 | 80.0 |
| mcr1 | Colistin | 42 | 41 (97.6%) | 1 | 86.9 | 84.8 | 78.8 |
| mecA | Methicillin/MRSA | 23 | 23 (100%) | 0 | 86.0 | 85.3 | 80.0 |
| **Total** | | **145** | **140 (96.6%)** | **5** | **87.8** | | |

**Key finding:** All top guides show zero off-target hits at 0–3 mismatches — 100% genome-wide specificity.

---

## Biological Motivation

Antimicrobial resistance kills approximately 1.27 million people annually and is projected to cause 10 million deaths per year by 2050 if unaddressed. The CRISPR-Cas9 system, reprogrammed with custom guide RNAs, can sequence-specifically disrupt resistance genes on both plasmids and bacterial chromosomes — re-sensitizing MDR bacteria to antibiotics that were previously ineffective.

Unlike broad-spectrum antibiotics, CRISPR-based antimicrobials offer programmable, gene-level precision, selectively eliminating resistance-carrying bacteria while sparing susceptible strains and commensal microbiota.

---

## Scoring Pipeline

### Composite Scoring Formula

```
Final Score = (On-Target Score × 0.6) + (Specificity Score × 0.4)

Specificity Score = max(0, 100 − Penalty)
Penalty = (0mm hits × 20) + (1mm hits × 5) + (2mm hits × 2) + (3mm hits × 0.5)
```

For pan-strain guides, conservation score is additionally incorporated:

```
Pan-Guide Final Score = (On-Target × 0.4) + (Specificity × 0.3) + (Conservation × 0.3)
```

### Classification Thresholds

| Classification | Final Score |
|---------------|-------------|
| Excellent | ≥ 85 |
| Good | 70 – 84 |
| Moderate | 50 – 69 |
| Poor | < 50 |

> Note: The pan-strain scoring scale ranges 76–88, reflecting the more stringent multi-strain filtering criteria.

---

## Pipeline Workflow

```
Gene Sequences (NCBI RefSeq)
        ↓
PAM Site Identification (NGG: TGG / CGG / AGG / GGG)
        ↓
20-nt Spacer Extraction
        ↓
On-Target Scoring (GC content 40–70%, Doench 2016 rules)
        ↓
Off-Target Analysis (0–3 mismatch genome-wide search)
        ↓
[Pan-Strain] Multi-Variant Conservation Scoring
        ↓
Composite Scoring & Classification
        ↓
Ranked Guide RNA Library
```

---

## Repository Structure

```
crispr-mdr-resensitization/
├── data/
│   └── targets_multistrain/        # Multi-strain FASTA inputs per gene
├── results/                        # Phase 1 single-reference outputs
│   ├── all_guide_candidates.csv
│   ├── top20_global_guides.csv
│   ├── top10_per_gene.csv
│   └── summary_statistics.csv
├── results_panstrain/              # Phase 2 pan-strain outputs
│   ├── all_panstrain_guide_candidates.csv
│   ├── top30_global_panstrain_guides.csv
│   ├── top20_per_gene_panstrain.csv
│   └── summary_statistics_panstrain.csv
├── figures/                        # All publication figures
│   ├── Figure1_ScoreDistribution.png
│   ├── Figure2_TopGuides_Global.png
│   ├── Figure3_GeneComparison.png
│   ├── Figure4_OnTarget_vs_OffTarget.png
│   ├── Figure1_PanStrainScoreDistribution.png
│   ├── Figure2_TopPanStrainGuides_Global.png
│   ├── Figure3_PanStrainGeneComparison.png
│   ├── Figure4_Specificity_vs_Conservation.png
│   ├── Figure5_GuideCoverageByGene.png
│   ├── Heatmap_blaKPC_TopGuides.png
│   ├── Heatmap_blaNDM1_TopGuides.png
│   ├── Heatmap_mcr1_TopGuides.png
│   └── Heatmap_mecA_TopGuides.png
├── crispr_mdr_analysis.py          # Main analysis script
├── .gitignore
└── README.md
```

---

## Output Column Descriptions

| Column | Description |
|--------|-------------|
| `gene` | Target resistance gene |
| `position` | Position in gene sequence |
| `strand` | + (sense) or − (antisense) |
| `spacer` | 20-nt guide RNA spacer sequence |
| `pam` | PAM sequence (TGG/CGG/AGG/GGG) |
| `gc_content` | GC% of spacer |
| `on_target_score` | On-target efficiency score (0–100) |
| `offtarget_hits_0mm` | Exact genome matches |
| `offtarget_hits_1mm` | 1-mismatch hits |
| `offtarget_hits_2mm` | 2-mismatch hits |
| `offtarget_hits_3mm` | 3-mismatch hits |
| `offtarget_penalty` | Weighted off-target penalty |
| `specificity_score` | 100 − penalty |
| `conservation_score` | Pan-strain conservation score (pan-strain files only) |
| `final_score` | Composite weighted score |
| `classification` | Excellent / Good / Moderate / Poor |

---

## Top Guides — Phase 1 (Score 90.4, 100% Specificity)

| Gene | Position | Spacer (5'→3') | PAM | GC% |
|------|----------|----------------|-----|-----|
| blaKPC | 756 | `CACAATAGGTGCGCGCCCAG` | TGG | 65% |
| blaNDM1 | 677 | `CGGTGATATTGTCACTGGTG` | TGG | 50% |
| mcr1 | 2650 | `AACCTGCGCAGGTAACACAG` | TGG | 55% |
| mecA | 1965 | `AAAGTGGCAGACAAATTGGG` | TGG | 45% |

## Top Guides — Phase 2 Pan-Strain (Score 87.8, 100% Specificity, 100% Strain Coverage)

| Gene | Position | Spacer (5'→3') | PAM | Conservation |
|------|----------|----------------|-----|-------------|
| blaNDM1 | 221 | `AATGGTTCGGTATGCGGGCG` | TGG | 80.0 |
| blaNDM1 | 308 | `AATGGTTCGGTATGCGGGCG` | TGG | 80.0 |
| blaKPC | 171 | `AGAGCCTTACTGCCCGAAGG` | CGG | 80.0 |
| blaKPC | 174 | `GCCTTACTGCCCGAAGGCGG` | CGG | 80.0 |

---

## Known Limitations & Planned Fixes

| Issue | Status |
|-------|--------|
| **GC content not penalized in Phase 1 on-target score** — 1,877 guides outside 40–70% range included in scored pool | Fix in progress — Phase 2 addresses this via pre-filtering |
| **mcr1 specificity score inconsistency (n=279)** — formula mismatch in Phase 1 processing | Under investigation |
| **No minimum on-target score threshold in Phase 1** — low-scoring guides inflate total candidate count | Planned: minimum threshold ≥ 50 |
| **Experimental validation pending** — all results are computational | Planned: CFPS in vitro validation, phage delivery testing |

---

## Requirements

```
python >= 3.8
biopython
pandas
numpy
matplotlib
seaborn
```

Install dependencies:
```bash
pip install -r requirements.txt
```

---

## Citation

If you use this pipeline or data, please cite:

> Basu, S. *CRISPR-Cas mediated re-sensitization of multidrug-resistant bacteria through plasmid and chromosomal targeting.* ACS Spring 2026, BIOT Division of Biochemical Technology, Paper ID: 4421933. University of Illinois Springfield.

---

## Contact

**Sanghati Basu**
University of Illinois Springfield
GitHub: [@SANGHATI23](https://github.com/SANGHATI23)
ACS Spring 2026 | Paper ID: 4421933

---

*This work is part of ongoing computational research into CRISPR-based antimicrobials for combating the global AMR crisis. All results are pre-publication and presented at ACS Spring 2026.*
