# BRCA-PRISM: Breast Cancer Pathway Resource Integration and Stability Mapping

BRCA-PRISM is a framework for evaluating how pathway-resource choice affects
gene-set composition, patient-level pathway activity, predictive performance,
and biological interpretation in breast cancer multi-omics analysis.

## Overview

This project compares five pathway/gene-set resources:

- KEGG
- Reactome
- WikiPathways
- MSigDB Hallmark
- LCPathways

The resources were harmonized into a common pathway-gene representation and
organized into cross-resource pathway families for downstream breast cancer analysis.

## Datasets

### 1. TCGA-BRCA
Primary breast cancer cohort used for pathway-resource comparison,
PAM50 subtype association, predictive analysis, and interpretation stability.

- Samples: 974 tumors
- Omics: RNA expression + CNV
- Labels: PAM50 molecular subtypes

### 2. CPTAC-BRCA
Independent breast cancer cohort used for external biological validation
and cross-cohort predictive sensitivity analysis.

- Samples: 122 tumors
- Omics: RNA expression + CNV
- Molecular subtype information available

### 3. TCGA-UCEC
Exploratory cross-disease cohort used to evaluate whether selected
BRCA-PRISM pathway-family signals extend beyond breast cancer.

- Samples: 507 tumors
- Cancer type: Uterine Corpus Endometrial Carcinoma
- Omics: RNA expression + CNV
- Molecular subtypes:
  - POLE-ultramutated
  - MSI-hypermutated
  - Copy-number-low
  - Copy-number-high
- Evaluated pathway families:
  - F0006 — DNA replication
  - F0025 — Endometrial cancer

## BRCA-PRISM Framework

The workflow includes:

1. Pathway database harmonization
2. Gene-level overlap analysis
3. Cross-resource pathway matching
4. Construction of BRCA-PRISM pathway families
5. TCGA-BRCA pathway activity analysis
6. PAM50 association analysis
7. Prediction-stability analysis
8. Interpretation-stability analysis
9. Consensus vs resource-specific gene ablation
10. External CPTAC-BRCA validation

## Key Results

- 231,507 pathway-gene associations
- 4,307 source pathways/gene sets
- 31 BRCA-PRISM cross-resource pathway families
- All 31 families were significantly associated with PAM50 subtype in TCGA-BRCA
- DNA replication and cell-cycle families showed the strongest subtype associations
- Gene-set similarity strongly predicted patient-level pathway-score concordance
- Similar predictions did not necessarily produce similar biological explanations
- Resource-specific genes performed better in TCGA cross-validation
- Consensus-family features showed stronger cross-cohort transfer behavior in the
  post-hoc TCGA-to-CPTAC sensitivity analysis
