🧬 CRISPR-MDR: Robust AMR Prediction to CRISPR Intervention Pipeline


 Overview
Antimicrobial resistance (AMR) is driven by diverse and rapidly evolving genetic mechanisms. While machine learning models can accurately predict resistance phenotypes from genomic data, and CRISPR-based systems enable programmable targeting of resistance genes, these two domains are typically treated independently.
This repository presents a unified computational framework that connects:
Genotype → Phenotype Prediction → Resistance Drivers → CRISPR Target Prioritization
Rather than stopping at prediction, this pipeline translates resistance-associated signals into actionable CRISPR intervention candidates, enabling computational prioritization of antimicrobial resensitization strategies.


 Core Contribution
This work introduces an integrated and robustness-tested pipeline that:
Predicts AMR phenotypes from genomic features
Identifies key resistance driver genes
Maps these genes to CRISPR targetable loci
Prioritizes guide RNAs using a multi-objective scoring framework
Demonstrates stability of guide selection across scoring models and weighting schemes


 Key Idea
Resistance prediction is not the endpoint — it becomes the input for intervention design.


 Conceptual Pipeline
Genomic data
   ↓
AMR feature matrix
   ↓
Machine learning model
   ↓
Resistance prediction
   ↓
Feature importance
   ↓
Key driver genes (mecA, blaKPC, blaNDM, mcr)
   ↓
CRISPR guide generation
   ↓
Multi-objective scoring
   ↓
Stable prioritized guides


 Methodology

1. AMR Prediction
Models: Logistic Regression, Random Forest, XGBoost
Evaluation: ROC-AUC, F1-score, Balanced Accuracy
Dataset aligned with CAMDA AMR challenge (official combinations)

2. Feature Importance → Biological Drivers
Top genes identified per pathogen-drug combination:
Gene	Role
mecA	Methicillin resistance (PBP2a)
blaKPC	Carbapenemase
blaNDM	Carbapenemase
mcr-1	Colistin resistance

3. Feature Importance → Biological Drivers
Top genes identified per pathogen-drug combination:
Gene	Role
mecA	Methicillin resistance (PBP2a)
blaKPC	Carbapenemase
blaNDM	Carbapenemase
mcr-1	Colistin resistance

4. Multi-objective Scoring
Final composite score:
Final Score = 0.40 × Activity + 0.30 × Specificity + 0.30 × Conservation
Where:
Activity → Doench/RS3-style efficiency
Specificity → CFD/MIT mismatch-aware penalty
Conservation → cross-strain robustness

5. Robustness Analysis (Key Strength)
66 weight configurations tested
Compared against:
On-target-only models
Off-target-only models
Conservation-based ranking
Consensus models
External comparators (CFD, MIT, CRISOT, RS3-style)

6. Weight Sensitivity Analysis (Robustness)
Instead of relying on a single scoring configuration, we evaluated 66 combinations of:
Activity (A)
Specificity (S)
Conservation (C)
Example:
40/30/30
30/30/40
30/40/30
50/30/20
...
Evaluation Metrics
Top-10 / Top-20 overlap
Jaccard similarity
Key Finding
Top-ranked guides remained stable across all configurations
Guide prioritization is not sensitive to weight tuning

7. Cross-Model Benchmarking
We compared the composite framework against multiple scoring paradigms:
| Category     | Models                  |
| ------------ | ----------------------- |
| On-target    | Doench 2016 / RS3-style |
| Off-target   | CFD, MIT                |
| Conservation | strain-aware proxy      |
| Ensemble     | consensus               |
| External     | CRISOT                  |

Results
Composite vs On-target → 10/10 overlap
Composite vs CFD/MIT → 8/10 overlap
Composite scoring balances efficiency + specificity + robustness


8. Multi-objective Scoring
Final composite score:
Final Score = 0.40 × Activity + 0.30 × Specificity + 0.30 × Conservation
Where:
Activity → Doench/RS3-style efficiency
Specificity → CFD/MIT mismatch-aware penalty
Conservation → cross-strain robustness

9. Composite Scoring Framework
Final scoring formula:
Final Score = 0.40 × Activity + 0.30 × Specificity + 0.30 × Conservation
Interpretation
Activity → cutting efficiency
Specificity → off-target risk
Conservation → strain robustness
Important
40/30/30 is not optimized, it is a balanced reference
Stability across weights is the key result

10. Robustness Analysis (Key Strength)
66 weight configurations tested
Compared against:
On-target-only models
Off-target-only models
Conservation-based ranking
Consensus models
External comparators (CFD, MIT, CRISOT, RS3-style)

 Key Results
🔹 ML Performance (CAMDA-aligned)
ROC-AUC up to 0.98+
Strong performance across pathogen-drug combinations

🔹 Guide Stability (Main Finding)
Top-ranked guides remained unchanged across all 66 weighting configurations
Top-20 overlap ≈ 100%
Cross-model agreement:
Composite vs On-target → 10/10 overlap
Composite vs CFD/MIT → 8/10 overlap

 Conclusion:
Guide prioritization is robust and model-invariant

🔹 CRISPR Targetability
Gene	Guides	Elite Guides
blaKPC	136	34
blaNDM-1	157	46
mcr-1	5962	42
mecA	118	23


 Repository Structure

results/                     → Core CRISPR scoring outputs  
results_external/           → CFD, MIT, CRISOT comparisons  
results_panstrain/          → Cross-strain conservation analysis  
results_panstrain_staph_200/→ Staphylococcus-specific validation  
results_weight_analysis/    → Weight sensitivity experiments  
weight_sensitivity_results/ → Stability metrics and overlaps  


 Figures (Manuscript)
AMR prediction performance (ROC curves)
Feature importance (driver genes)
CRISPR score distributions
Weight sensitivity robustness
Cross-model comparisons (CFD, MIT)

 Integrated pipeline (Prediction → Intervention bridge)

 Limitations

Fully computational (no wet-lab validation)
Gene presence ≠ expression or phenotypic dominance
Conservation approximated via available strain data

 Future Directions

Experimental validation of top-ranked guides
Integration with expression data (transcriptomics)
Expansion to broader pathogen-drug space
CRISPR delivery strategy modeling

 Final Positioning
This work demonstrates that CRISPR guide prioritization for AMR targets can be made robust across scoring assumptions, and establishes a direct computational bridge from genomic prediction to intervention design.

 Citation (if used)
Basu, S. (2026)
A Robust and Interpretable Framework Linking AMR Prediction to CRISPR-Based Intervention Prioritization

 Contact
For collaboration or discussion:
Sanghati Basu
University of Illinois Springfield

 FINAL NOTE (THIS IS WHAT PIs SEE)
This is not just:
ML for AMR
OR CRISPR scoring
 It is:
A unified, robustness-validated framework that connects prediction to intervention.
