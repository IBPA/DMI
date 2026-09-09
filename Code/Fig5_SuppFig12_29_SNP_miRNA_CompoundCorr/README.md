# SNP–compound and miRNA–compound correlation analyses used for Figure 5 and Supplementary Figures 12 and 29

SNP/miRNA-to-compound correlation analyses and supporting scripts.

## Scripts

### `SNP_CompoundCorrelation.r`

Summarizes significant SNP–compound correlations.

**Input:** Either of the following:

1. An SNP–compound correlation matrix with corresponding p-values. These large files are available upon request.
2. Prepared significant SNP–compound pairs (example: `ProcessedData/Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv`).

The variable `USE_PREPARED_SNP_COMPOUND_PAIRS` determines whether to use the prepared pairs as input.

The script also requires bioactive compound data with estimated bioactivity efficacy, comprising the merged estimates of in vitro and in vivo bioactivity efficacy of bioassays in milk samples. Files are available upon request.

**Output:** A summary of SNP–compound pairs for Supplementary Figure 29a and Figure 5. Additional manual preparation is required to produce the publication figures.

### `miRNA_CompoundCorrelation.r`

Summarizes miRNA–compound correlations.

**Input:**

1. Concentration data and miRNA counts prepared by `Code/1_ConcentrationTablePreparation.R`.
2. Bioactive compound data with estimated bioactivity efficacy, comprising the merged estimates of in vitro and in vivo bioactivity efficacy of bioassays in milk samples. Files are available upon request.

**Output:** A summary of miRNA–compound pairs for Supplementary Figures 12 and 29b. Additional manual preparation is required to produce the publication figures.

## Subdirectories

- `Python_CorrelationEvaluation/`: Contains Python code for calculating SNP–compound correlation matrices using Pearson correlation coefficients (PCC) and Spearman correlations, along with corresponding p-values, from SNP allele frequencies and compound concentration data. These input data are available upon request.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
