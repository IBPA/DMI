# Bioactivity analyses used for Figure 3 and Supplementary Figures 6 and 17

## Scripts

### `GetSmallMoleculeBioactivity.R`

Retrieves PubChem bioassays that tested small molecules identified in milk.

**Input:** Concentration data prepared by `Code/1_ConcentrationTablePreparation.R.

Internet access is required to download PubChem bioassay data.

**Output:**

1. Bioassays with Hill curves. (Example: `ProcessedData/Fig3_SuppFig6_17_Bioactivity/MilkSmallMoleculeBioactivity_FittedCurve.csv`.)
2. Descriptions of selected PubChem bioassays related to health. Additional manual labeling is performed for the next step. (Examples: `ProcessedData/Fig3_SuppFig6_17_Bioactivity/BioassayDescription.csv` and `ProcessedData/Fig3_SuppFig6_17_Bioactivity/BioassayDescription_Selected_Manual.csv`.)

### `SmallMoleculeBioactivityEfficacyEstimation.R`

Estimates both in vitro and in vivo bioactivity of the milk samples.

**Input:**

1. Bioassays with Hill curves. (Example: `ProcessedData/Fig3_SuppFig6_17_Bioactivity/MilkSmallMoleculeBioactivity_FittedCurve.csv`.)
2. Descriptions of selected PubChem bioassays related to health, with additional manual labeling. (Examples: `ProcessedData/Fig3_SuppFig6_17_Bioactivity/BioassayDescription.csv` and `ProcessedData/Fig3_SuppFig6_17_Bioactivity/BioassayDescription_Selected_Manual.csv`.)

**Output:** Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to small molecules with Hill curves. Files are available upon request.

### `PeptideEfficacyEstimation.R`

Estimates in vitro and in vivo bioactivity of the milk samples attributable to peptides.

**Input:**

1. Manually curated bioactive peptide data with IC50 values. (Example: `Data/Bioactivity/ManuallyCuratedPeptideIC50Data.csv`.)
2. Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to small molecules with Hill curves. Files are available upon request.

**Output:** Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to peptides. Files are available upon request.

### `HillCurvePerformanceTest.R`

Evaluates the error introduced by using representative Hill curves to estimate bioactivity efficacy when only IC50 values are available.

**Input:** Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to small molecules with Hill curves. Files are available upon request.

**Output:** Bioactivity efficacy estimation error assessed by bootstrapping (Supplementary Figure 17).

### `MergeSmallMoleculeAndPeptideBioactivity.R`

Merges the estimated in vitro and in vivo bioactivity efficacy attributable to small molecules and peptides.

**Input:**

1. Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to small molecules. Files are available upon request.
2. Estimated in vitro and in vivo bioactivity efficacy for selected bioassays in milk samples, attributable to peptides. Files are available upon request.

**Output:** Merged estimated in vitro and in vivo bioactivity efficacy of bioassays in milk samples. Files are available upon request.

### `PlotBioactivityEfficacy.R`

Plots the bioactivity efficacy estimation results.

**Input:** Merged estimated in vitro and in vivo bioactivity efficacy of bioassays in milk samples. Files are available upon request.

**Output:** Figure 3 and Supplementary Figure 6.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
