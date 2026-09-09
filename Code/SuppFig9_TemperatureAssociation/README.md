# Temperature and packaging association analyses used for Supplementary Figure 9

## Script

### `Temperature_Packaging_Interaction.R`

Examines associations between selected compound concentrations, ambient temperature at the purchase location, and packaging opacity, including temperature–packaging interactions.

**Input:**

- Concentration data prepared by `Code/1_ConcentrationTablePreparation.R`.
- Product information containing nutritional facts in `Data/Products.csv`.

**Output:** Cases in which ambient temperature at the purchase location and packaging opacity are associated with selected compound concentrations (Supplementary Figure 9).

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
