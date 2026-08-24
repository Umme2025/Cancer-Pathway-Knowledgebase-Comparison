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

## Main Dataset

- TCGA-BRCA: 974 tumors
- Molecular data: RNA expression + CNV
- Outcome: PAM50 breast cancer subtype
- Independent validation: CPTAC-BRCA, n = 122

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
