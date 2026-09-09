# Code Directory

This directory contains analysis workflows organized by figure and supplemental figure group.

## Contents

- `Fig2a_SuppFig2_3_4_16_CompoundComposition/`: Compound composition analyses used for Figure 2a and related supplemental figures.
- `Fig2bc_SuppFig5_NutritionalFactComparison/`: Nutritional fact comparison analyses used for Figure 2b/c and Supplemental Figure 5.
- `Fig3_SuppFig6_17_Bioactivity/`: Bioactivity data integration, curation, efficacy estimation, and plotting workflows.
- `Fig4_SuppFig8_24_RegionAnalysis/`: Region-level analyses including map plotting and hierarchical variance decomposition.
- `Fig5_SuppFig12_29_SNP_miRNA_CompoundCorr/`: SNP/miRNA-to-compound correlation analyses and supporting scripts.
- `NNLS-SNP-main/`: Packaged NNLS-SNP breed/variety composition estimator module and source code. The reference allele frequency file curated from the BGVD database is available upon request.
- `SuppFig10_11_JerseyBreedAssociation/`: Jersey breed association analyses, ratio estimation, and supporting plots.
- `SuppFig14_SNPSummary/`: SNP summary analyses used for Supplemental Figure 14.
- `SuppFig15_CompoundPCA/`: Principal component analysis workflows for compound concentration data.
- `SuppFig18To23_MOFA/`: MOFA (multi-omics factor analysis) workflows for Supplemental Figures 18–23.
- `SuppFig1_2_SuppTable1_OmicsSetComparison/`: Comparison workflows across omics datasets and external databases.
- `SuppFig25To28_SNPConstructCompoundProfiles/`: SNP-based compound profile construction analyses for Supplemental Figures 25–28.
- `SuppFig30_MediationAnalysis/`: Mediation analysis workflows and plotting scripts for Supplemental Figure 30.
- `SuppFig9_TemperatureAssociation/`: Temperature and packaging interaction analysis workflows.

## Additional notes

- `NNLS-SNP-main/` documentation is available in the dedicated repository: https://github.com/IBPA/NNLS-SNP
- Python code is present in:
  - `Fig5_SuppFig12_29_SNP_miRNA_CompoundCorr/Python_CorrelationEvaluation/` (for generating PCC/Spearman correlation matrices)
  - `SuppFig15_CompoundPCA/` (entire directory is Python-based)
- For these Python workflows, necessary input data files are stored in the same directory as the Python code and can be provided upon request.
