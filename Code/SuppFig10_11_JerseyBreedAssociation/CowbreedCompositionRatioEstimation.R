#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

dir.create(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation"), showWarnings = FALSE)

source("Code/NNLS-SNP-main/NNLS-SNP-main/src/nnls_functions.R")
library(stringr)

load(paste0(DATA_PATH,"/8MSNPs.rdata"))
load(paste0(CODE_PATH,"/NNLS-SNP-main/NNLS-SNP-main/data/three_breeds_allele_freq.rdata"))
ad_freq = ad_only_one_alt_2_intersect/ad_only_one_alt_sum_intersect

#ref_allele_frequency
new_rownames = sapply(strsplit(rownames(ad_freq),"_"), function(x){
  A = gsub("GK0000","",x[1])
  A = as.numeric(gsub("[.]2","",A))
  if (is.na(A)) return(NA)
  
  return(paste0(A,":",x[2],"_",x[3],"/",x[4]))
})

rownames(ad_freq) = new_rownames


val_data_three_breeds_matrix_numeric = apply(val_data_three_breeds_matrix, 2, as.numeric)
rownames(val_data_three_breeds_matrix_numeric) = rownames(val_data_three_breeds_matrix)
intersect_snps = intersect(rownames(ad_freq), rownames(val_data_three_breeds_matrix_numeric))


result = run_nnls_merge_batch(as.matrix(ad_freq[intersect_snps,]),
                              as.matrix(val_data_three_breeds_matrix_numeric[intersect_snps,c(2,3)]))
write.csv(result,paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv"))
save(ad_freq, val_data_three_breeds_matrix_numeric, intersect_snps, file=paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_adfreq_matrix.rdata"))
