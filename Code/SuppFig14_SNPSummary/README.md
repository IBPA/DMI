# SNP summary analyses used for Supplementary Figure 14

## Script

### `SNP_Summary.R`

Summarizes SNP impact levels, detection frequency across samples, and sequencing depth.

**Input:**

1. An allele frequency matrix for 8 million SNPs. The file is available upon request.
2. An SNP outcome-to-impact-level mapping table (example: `Data/SNPInformation/SNPEffImpactLevel.csv`).

**Output:**

- **Supplementary Figure 14b:** Numbers of SNPs at different impact levels.
- **Supplementary Figure 14c:** Proportions of SNPs detected in fewer than 90%, 90% to less than 95%, and at least 95% of samples.
- **Supplementary Figure 14d:** Proportions of SNPs with a sequencing depth of at least 5 in fewer than 50%, 50% to less than 80%, and at least 80% of samples.

**Manual preparation of Supplementary Figure 14a:** This panel was prepared by manually comparing detected SNPs with SNPs recorded in the BGVD database. The allele frequency matrix for all detected SNPs (approximately 26 million) is available upon request.

## Usage

Run scripts in this folder as part of the corresponding figure workflow, using shared paths configured in `Code/FilePath.R` where applicable.
