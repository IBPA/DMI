# Multi-omics factor analysis used for Supplementary Figures 18–23

## Script

### `RunMoFA.R`

Performs multi-omics factor analysis (MOFA) using compound concentration data, miRNA counts, and SNP allele frequencies.

**Input:** The following data files are available upon request:

1. Concentration data and miRNA counts prepared by `Code/1_ConcentrationTablePreparation.R`.
2. An allele frequency matrix for SNPs with medium impact levels.

**Output:**

- **Supplementary Figure 18:** Scatter plots of MOFA factor scores for Factor 1 versus Factor 2 and Factor 2 versus Factor 3.
- **Supplementary Figure 19:** A summary of the MOFA analysis, including explained variance and absolute scaled weights for the first two factors.
- **Supplementary Figures 20 and 21:** Absolute weights of selected features for Factor 1 in each omics layer.
- **Supplementary Figures 22 and 23:** Absolute weights of selected features for Factor 2 in each omics layer.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
