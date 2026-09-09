# SNP-based compound profile reconstruction used for Supplementary Figures 25–28

## Script

### `SNPConstructCompoundProfiles.R`

Predicts principal component (PC) scores of compound profiles from selected SNPs and evaluates compound profile reconstruction.

**Input:**

1. Concentration data and miRNA counts prepared by `Code/1_ConcentrationTablePreparation.R`.
2. An allele frequency matrix for SNPs with medium impact levels.

**Output:**

- **Supplementary Figure 25:** Scatter plots comparing actual compound profile PC scores with scores predicted from 61 selected SNPs. The linear models are trained on the full dataset.
- **Supplementary Figure 26:** Linear model coefficients for predicting seven PCs, with the models trained on the full dataset.
- **Supplementary Figure 27:** Leave-one-region-out cross-validation performance for PC score prediction and compound profile reconstruction.
- **Supplementary Figure 28:** Compound profile reconstruction performance, measured by Pearson correlation coefficients (PCC) and Spearman correlations, evaluated by bootstrapping with different proportions of test data.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
