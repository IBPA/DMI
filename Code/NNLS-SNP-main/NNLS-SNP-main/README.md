
# Breed/Variety Composition Ratio Estimator
This module evaluates the breed or variety composition ratio of samples based on the SNP (Single Nucleotide Polymorphism) allele frequency information.

![Figure 1. The overview of the sequence novelty evaluator.](https://github.com/IBPA/NNLS-SNP/blob/main/Flow.png)

Two inputs must be provided to evaluate the composition ratio:

 - The SNP allele frequency of the samples.
	 - Must be a matrix (in R environment) with specified `rownames` and `colnames`, or in a VCF file format that stores the allele frequency of ALL samples (See **Section 2** for details).
 - The SNP allele frequency of the reference breed or variety (Curated from the published database, such as the [BGVD database](http://animal.omics.pro/code/index.php/BosVar) for cow breed allele frequency data)
	 - Must be a matrix (in R environment) with specified `rownames` and `colnames` (See **Section 2** for details).

## 1. Directory
There is a simple sub-repository. All input files are stored in the `data` directory.

 - Input:
	 - `cowbreed_samplemix_allele_frequency.rdata`:  The SNP allele frequency of 16 synthetic cow breed mixture samples.
		 - Stored as a variable `samplemix_allele_frequency`
		 - Hereford, Jersey, and Holstein cow breeds are mixed in different ratios.
		 - The sample label (`colnames`) recorded the mix ratios, for example:
			 - The Label `H10_J20_Hol5` means the ratio of Hereford, Jersey, and Holstein is (10 : 20 :5)
			 - The Label `H10_JX_Hol5` means the ratio of Hereford, Jersey, and Holstein is (10 : 0 : 5)
	 - `cowbreed_ref_allele_frequency.rdata`: The SNP allele frequency of the reference cow breeds (Hereford, Jersey, and Holstein in this case)
		 -  Stored as a variable `ref_allele_frequency`

## 2. Getting Started
The code is tested under the R environment version 4.2.0 with the following packages: `nnls`, `foreach`, and `doParallel`. You may install them using the following command in R:

```r (4.2.0)
install.packages(c("nnls", "foreach", "doParallel"))
```
The source code file `nnls_functions.R` contains several functions. To use that, just source the file:
```r (4.2.0)
source("src/nnls_functions.R")
```

### (Optional) Preparing the Reference Allele Frequency Data
The reference allele frequency data is a matrix with specified `rownames` and `colnames`.
 - `rownames`: The SNP name or ID. (**Note:** The SNP name or ID should be consistent with the SNP name or ID in the sample allele frequency data)
 - `rownames`: The breed (variety) names.

The allele frequency values are between 0 to 1. Here is an example (only showing the first six SNPs):
```r
> head(ref_allele_frequency)
            Hereford_val Holstein_val Jersey_val
1:324_A/G          0.167        0.300      0.273
1:340_G/A          0.167        0.300      0.250
1:380_G/T          0.167        0.300      0.208
1:2377_A/G         0.048        0.311      0.375
1:11412_C/G        0.476        0.433      0.083
1:14259_C/T        0.500        0.500      0.458
```
### (Optional) Preparing the Allele Frequency Data of the Mixture Samples
The allele frequency data of the mixture samples is a matrix with specified `rownames` and `colnames`.
 - `rownames`: The SNP name or ID. (**Note:** The SNP name or ID should be consistent with the SNP name or ID in the reference allele frequency data)
 - `rownames`: The sample names.

The allele frequency values are between 0 to 1. Here is an example (only showing the first six SNPs and the first five samples):
```r
> head(samplemix_allele_frequency)
            H10_J20_Hol5 H10_J5_Hol20 H10_J5_HolX H10_JX_Hol5 H20_J10_Hol5
1:324_A/G      0.7941176    0.2000000   0.6363636   0.0000000    0.4210526
1:340_G/A      0.9090909    0.9166667   0.7272727   0.5000000    0.4324324
1:380_G/T      0.7857143    0.8604651   0.6000000   0.3333333    0.3888889
1:2377_A/G     1.0000000    1.0000000   1.0000000   1.0000000    1.0000000
1:11412_C/G    0.1481481    0.1666667   0.5714286   0.4000000    0.4000000
1:14259_C/T    0.3157895    0.2500000   0.6666667   1.0000000    0.6285714
```
#### Tools for preparing the Allele Frequency Data Matrix from the VCF file
The `extract_ad_matrix_info_batch` (in `nnls_functions.R`) can help you to get the Allele Frequency Data Matrix from the VCF file generated from the tools for SNP detection (e.g., [GATK](https://gatk.broadinstitute.org/hc/en-us)) 

(**Note**: The VCF files must contain the allele frequencies of ALL samples, so they should be yielded from the joint-genotyping approach or combined from multiple VCF files that contain the allele frequencies of a single sample):

```r
extract_result = extract_ad_matrix_info_batch("Batch1-gatk-haplotype-annotated.vcf.gz")
samplemix_allele_frequency = extract_result$ad_matrix
```

### 2. Composition Ratio Evaluation:

Run the function `run_nnls_merge_batch()` to evaluate the breed/variety composition ratio of the samples:
```r
#Load the input data
>load('cowbreed_ref_allele_frequency.rdata')
>load('cowbreed_samplemix_allele_frequency.rdata')
>result = run_nnls_merge_batch(samplemix_allele_frequency, ref_allele_frequency)
>head(result)
             Hereford_val Holstein_val Jersey_val
H10_J20_Hol5    0.1992635   0.22503754  0.5756990
H10_J5_Hol20    0.1817231   0.69454493  0.1237320
H10_J5_HolX     0.6538215   0.04202237  0.3041561
H10_JX_Hol5     0.6180149   0.38198506  0.0000000
H20_J10_Hol5    0.5362635   0.20531296  0.2584236
H20_J5_Hol10    0.5307980   0.36143639  0.1077656
```

## 3. Authors

 - **ChengEn Tan**  @https://github.com/bigghost2054

## 4. Contact

For any questions, please contact us at  [tagkopouloslab@ucdavis.edu](mailto:tagkopouloslab@ucdavis.edu).

## 5. Citation

 - *(Will be ready once the work is published).*

## 6. License

The entire project is licensed under the Apache 2.0 License. Please see the  [LICENSE](https://github.com/IBPA/NNLS-SNP/blob/main/LICENSE)  file for details.


> Written with [StackEdit](https://stackedit.io/).