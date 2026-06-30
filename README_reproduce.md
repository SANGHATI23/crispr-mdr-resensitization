# FOTR-CRISPR Reproducibility Guide

This repository contains the reproducible computational workflow for the FOTR-CRISPR framework: a functional-context and RNA-structure-aware guide prioritization approach for antimicrobial resistance genes.

## Project scope

FOTR-CRISPR prioritizes CRISPR-Cas9 guide candidates against antimicrobial resistance genes by integrating:

- guide activity
- specificity
- pan-strain conservation
- full sgRNA RNA structural accessibility
- target-side functional context
- risk penalties
- external CRISOT spacer-level comparison
- ablation analysis
- mecA case-study interpretation

This workflow is intended for computational prioritization only. It does not claim experimental cleavage, antimicrobial resensitization, clinical efficacy, or delivery validation.

## Environment

The workflow was developed in VS Code on macOS using a Python virtual environment.

Required external dependency:

```bash
RNAfold