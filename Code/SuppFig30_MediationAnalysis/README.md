# Mediation analysis used for Supplementary Figure 30

## Script

### `MediationAnalysis.R`

Examines whether Jersey breed proportion mediates associations between SNPs and compounds.

**Input:**

1. Sample code mapping data in `Data/SampleCodeMapping.csv`.
2. Output files from `CowbreedCompositionRatioEstimation.R`, particularly the estimated Jersey and Holstein breed proportions for each sample.
3. Compounds correlated with Jersey breed proportion (example: `ProcessedData/SuppFig10_11_JerseyBreedAssociation/Compound_Correlated_to_Jersey_Breed.csv`).
4. Manually prepared mediation analysis results (example: `ProcessedData/SuppFig30_MediationAnalysis/lavaan_analysis_result_added_snp_metadata_manual.csv`), or all three of the following inputs: an allele frequency matrix for SNPs with medium impact levels, significant SNP–compound pairs, and bioactivity efficacy tables.

**Output:** Mediation analysis results for Supplementary Figure 30, including the numbers of SNP–compound pairs classified as mediated, not mediated, or partially mediated by Jersey breed proportion.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
