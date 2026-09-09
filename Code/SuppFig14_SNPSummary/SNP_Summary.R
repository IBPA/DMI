#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

load(paste0(DATA_PATH,"/8MSNPs.rdata"))
rownames(bgvd_data_intersect) = paste0(bgvd_data_intersect$UMD3_1_1_Pos,"_",bgvd_data_intersect$Alleles)

#Supp Fig 14 b
SNP_impact_level = read.csv(paste0(DATA_PATH,"/SNPInformation/SNPEffImpactLevel.csv"), header = T, check.names = F)
impact_level_dict = lapply(SNP_impact_level$Impact, function(x){x})
names(impact_level_dict) = SNP_impact_level$`Effect Seq. Ontology`
bgvd_data_intersect$impact_levels = sapply(bgvd_data_intersect$ConsequenceType, function(x){impact_level_dict[[x]]})
bgvd_data_intersect$impact_levels = stringr::str_to_title(bgvd_data_intersect$impact_levels)

ad_only_one_alt_sum_intersect_for_count = ad_only_one_alt_sum_intersect[,-c(7,8)]

detected_sample_percentage = apply(ad_only_one_alt_sum_intersect_for_count, 1, function(x){
  return(length(which(!is.na(x)))/ncol(ad_only_one_alt_sum_intersect_for_count))
})

#Supp Fig 14 c
detected_sample_lt_90 = length(which(detected_sample_percentage < 0.9))
detected_sample_ge_95 = length(which(detected_sample_percentage >= 0.95))
detected_sample_90_95 = length(which(detected_sample_percentage >= 0.90 & detected_sample_percentage < 0.95))

depth_count = apply(ad_only_one_alt_sum_intersect_for_count, 1, function(x){
  return(length(which(x >= 5))/ncol(ad_only_one_alt_sum_intersect_for_count))
})

#Supp Fig 14 d
depth_count_gt_80 = length(which(depth_count >= 0.8))
depth_count_50_80 = length(which(depth_count >= 0.5 & depth_count < 0.8))
depth_count_lt_50 = length(which(depth_count < 0.5))

#The pie chart in Sup Fig 14c and d is created based on the above statistics manually.