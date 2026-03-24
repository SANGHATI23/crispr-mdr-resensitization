# CRISPR-MDR Re-Sensitization

Computational CRISPR guide RNA design pipeline for antimicrobial resistance gene targeting in multidrug-resistant bacteria.

**ACS Spring 2026**  
**Division:** BIOT – Division of Biochemical Technology  
**Session:** Microbial Engineering and Fermentation  
**Paper ID:** 4421933

---

## Overview

This repository presents a computational pipeline for designing, scoring, ranking, and comparing CRISPR guide RNA candidates targeting major antimicrobial resistance (AMR) genes in multidrug-resistant bacteria.

The project focuses on four clinically important resistance genes:

- **blaKPC** – carbapenem resistance
- **blaNDM1** – carbapenem resistance
- **mcr1** – colistin resistance
- **mecA** – methicillin resistance / MRSA

The goal is to identify high-confidence CRISPR guide RNAs with strong predicted activity, low off-target risk, and improved pan-strain relevance for downstream experimental validation.

---

## Dataset Summary

A total of **6,373 guide RNA candidates** were screened across four major multidrug-resistance genes.

| Gene | Resistance Type | Total Guides | Excellent Guides | Best Score |
|------|------------------|-------------:|-----------------:|-----------:|
| blaKPC | Carbapenem | 136 | 129 (94.9%) | 90.4 |
| blaNDM1 | Carbapenem | 157 | 153 (97.5%) | 90.4 |
| mcr1 | Colistin | 5,962 | 5,004 (83.9%) | 90.4 |
| mecA | Methicillin / MRSA | 118 | 74 (62.7%) | 90.4 |

---

## Biological Motivation

Antimicrobial resistance is one of the most urgent challenges in modern infectious disease research. This project explores a computational CRISPR-based strategy for **re-sensitization**, where resistance-associated genes are targeted to potentially reduce resistance burden and support future therapeutic or experimental interventions.

The selected genes represent major resistance mechanisms in clinically important bacterial pathogens:

- **blaKPC**: carbapenemase-associated resistance
- **blaNDM1**: metallo-beta-lactamase-mediated resistance
- **mcr1**: plasmid-mediated colistin resistance
- **mecA**: methicillin resistance in staphylococcal strains

---

## Project Objectives

- Identify CRISPR-Cas9 guide RNA candidates against key AMR genes
- Rank guides using on-target, off-target, and conservation-aware scoring
- Compare top candidates across multiple resistance genes
- Extend analysis to pan-strain / multi-strain sequence contexts
- Generate publication-ready figures and summary tables

---

## Main Components

### `crispr_mdr_analysis.py`
Primary analysis script responsible for:
- loading input target sequences
- identifying PAM-compatible CRISPR guide candidates
- scoring and ranking guides
- generating summary tables
- producing visual outputs

### `data/targets_multistrain/`
Contains the multi-strain sequence inputs used for pan-strain analysis of:
- blaKPC
- blaNDM1
- mcr1
- mecA

### `results_panstrain/`
Contains pan-strain outputs, including:
- `all_panstrain_guide_candidates.csv`
- `top20_per_gene_panstrain.csv`
- `top30_global_panstrain_guides.csv`
- `summary_statistics_panstrain.csv`
- pan-strain heatmaps
- conservation and specificity plots
- ranked comparison figures

---

## Pipeline Workflow

The computational workflow follows these steps:

1. Load antimicrobial resistance gene sequences
2. Detect PAM-compatible CRISPR target sites
3. Extract candidate guide RNA sequences
4. Estimate predicted on-target performance
5. Evaluate off-target burden
6. Incorporate multi-strain conservation evidence
7. Compute composite guide scores
8. Rank and export top candidates
9. Generate visual summaries and comparison figures

---

## Scoring Strategy

Guide candidates are prioritized using a composite scoring logic that considers:

- **On-target efficiency**
- **Off-target specificity**
- **Pan-strain conservation / relevance**

This helps prioritize guide RNAs that are not only effective in a single sequence context, but also more robust across multiple strain variants.

---

## Classification

| Classification | Score |
|----------------|------:|
| Excellent | ≥ 80 |
| Good | 65 – 79 |
| Moderate | 50 – 64 |
| Poor | < 50 |

---

## Example Outputs

This repository generates outputs including:

- ranked guide candidate tables
- per-gene top guide summaries
- global top guide rankings
- score distribution plots
- pan-strain heatmaps
- conservation versus specificity figures
- comparative gene-level visualizations

These outputs support interpretation, prioritization, and presentation of computational findings.

---

## Repository Structure

```text
crispr-mdr-resensitization/
├── data/
│   └── targets_multistrain/
├── results/
├── results_panstrain/
├── .gitignore
├── crispr_mdr_analysis.py
├── guide_output.txt
└── README.md
