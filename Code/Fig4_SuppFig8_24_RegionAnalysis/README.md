# Regional analyses used for Figure 4 and Supplementary Figures 8 and 24

## Scripts

### `HierarchicalVarianceDecomposition.R`

Estimates the proportion of variance attributable to each factor (Region × Brand, Region, and Time) and the residual. It also evaluates regional effects on the concentrations of selected bioactive compounds by comparing estimated marginal means (EMMs) across milk samples from different regions.

**Input:**

1. Concentration data prepared by `Code/1_ConcentrationTablePreparation.R`.
2. Bioactive compound data with estimated bioactivity efficacy, comprising the merged estimates of in vitro and in vivo bioactivity efficacy of bioassays in milk samples. Files are available upon request.

**Output:**

1. Estimated proportions of variance attributable to Region × Brand, Region, Time, and the residual (Supplementary Figures 8 and 24).
2. Bioactive compounds whose concentrations differ significantly across regions, with comparisons of their concentrations and bioactivity efficacy (Figure 4).

### `Plot_USMap_IRIRegion.r`

Plots a U.S. map colored according to the regions defined by IRI. Additional manual preparation is required to produce the publication figures.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
