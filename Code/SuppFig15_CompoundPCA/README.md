# Compound principal component analysis used for Supplementary Figure 15

Principal component analysis workflows for compound concentration data.

## Python scripts

### `PCA.py`

Performs principal component analysis (PCA) of compound concentration profiles in milk samples.

**Input:** The following data files are available upon request:

1. A milk compound concentration matrix in CSV format.
2. Milk molecule metadata mapping compound IDs to compound names.
3. Milk product information recording the region code, brand code, and purchase time order for each sample.

**Output:** PCA results for the milk compound profiles (Supplementary Figure 15a).

### `pca_loading.py`

Identifies the compounds with the largest absolute loadings on the first two principal components (PCs).

**Input:** The same concentration matrix, molecule metadata, and product information used by `PCA.py`. These files are available upon request.

**Output:** Compounds with the largest absolute loadings on PC1 and PC2 (Supplementary Figure 15b).
