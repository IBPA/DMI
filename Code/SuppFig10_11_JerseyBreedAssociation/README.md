# Jersey breed association analyses used for Supplementary Figures 10 and 11

## Scripts

### `CowbreedCompositionRatioEstimation.R`

Estimates the proportions of Jersey and Holstein breeds in each milk sample.

**Input:**

1. An allele frequency matrix for 8 million SNPs. The file is available upon request.
2. Reference allele frequencies for Jersey and Holstein breeds collected from the BGVD database. The file can be prepared from BGVD data and is also available upon request.

**Output:** Estimated Jersey and Holstein breed proportions for each sample (example: `ProcessedData/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv`).

### `Heatmap_SNPs.R`

Plots heatmaps of SNP allele frequencies.

**Input:**

1. Output files from `CowbreedCompositionRatioEstimation.R`, particularly the estimated Jersey and Holstein breed proportions for each sample.
2. An allele frequency matrix for 8 million SNPs. The file is available upon request.

**Output:** Heatmaps of SNP allele frequencies (Supplementary Figure 10).

### `BreedRatioBoxplot.R`

Compares estimated Jersey breed proportions across regions and examines their association with kappa-casein concentration.

**Input:**

1. Output files from `CowbreedCompositionRatioEstimation.R`, particularly the estimated Jersey and Holstein breed proportions for each sample.
2. Sample code mapping data in `Data/SampleCodeMapping.csv`.
3. Concentration data prepared by `Code/1_ConcentrationTablePreparation.R`.

**Output:**

1. A boxplot of estimated Jersey breed proportions in milk samples from 10 regions.
2. The correlation between estimated Jersey breed proportion and kappa-casein concentration.

These results are used for Supplementary Figure 11.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
