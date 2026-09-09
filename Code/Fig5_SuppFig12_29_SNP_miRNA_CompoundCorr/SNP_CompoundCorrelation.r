#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))
library(stringr)

USE_PREPARED_SNP_COMPOUND_PAIRS = T

if (USE_PREPARED_SNP_COMPOUND_PAIRS){
  # Load the prepared significant snp-compound pairs (q-value threshold: 0.1; correlation threshold: 0.3; non-zero count threshold: 10)
  if (!file.exists(paste0(PROCESSED_DATA_PATH, "/Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv"))){
    stop("Missing prepared significant SNP-compound pairs data (Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv).")
  }
  
  significant_snp_compound_df = read.csv(paste0(PROCESSED_DATA_PATH, "/Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv"))
  
}else{
  if (!file.exists(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\pcc_correlation_matrix.csv.gz")) ||
      !file.exists(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\spearman_correlation_matrix.csv.gz")) ||
      !file.exists(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\pcc_pvalues.csv.gz")) || 
      !file.exists(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\spearman_pvalues.csv.gz"))){
    stop("Missing required Pearson and Spearman correlation matrix and p-value data.")
  }
  
  thres_corr = 0.3
  non_zero_count_thres = 10
  thres_qval = 0.1
  
  #Load Correlation Data
  pcc_snp_compound_correlation = read.csv(gzfile(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\pcc_correlation_matrix.csv.gz"), 'rt'), row.names=1, check.names = F)
  spearman_snp_compound_correlation = read.csv(gzfile(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\spearman_correlation_matrix.csv.gz"), 'rt'), row.names=1, check.names = F)
  pval_pcc_snp_compound_correlation = read.csv(gzfile(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\pcc_pvalues.csv.gz"), 'rt'), row.names=1, check.names = F)
  pval_spearman_snp_compound_correlation = read.csv(gzfile(paste0(DATA_PATH,"\\SNP_CompoundCorrelation\\spearman_pvalues.csv.gz"), 'rt'), row.names=1, check.names = F)
  
  qval_pcc_snp_compound_correlation = matrix(p.adjust(as.numeric(as.matrix(pval_pcc_snp_compound_correlation)), method = "BH"),
                                             nrow = nrow(pval_pcc_snp_compound_correlation),
                                             ncol = ncol(pval_pcc_snp_compound_correlation))
  rownames(qval_pcc_snp_compound_correlation) = rownames(pval_pcc_snp_compound_correlation)
  colnames(qval_pcc_snp_compound_correlation) = colnames(pval_pcc_snp_compound_correlation)
  
  qval_spearman_snp_compound_correlation = matrix(p.adjust(as.numeric(as.matrix(pval_spearman_snp_compound_correlation)), method = "BH"),
                                                  nrow = nrow(pval_spearman_snp_compound_correlation),
                                                  ncol = ncol(pval_spearman_snp_compound_correlation))
  rownames(qval_spearman_snp_compound_correlation) = rownames(pval_spearman_snp_compound_correlation)
  colnames(qval_spearman_snp_compound_correlation) = colnames(pval_spearman_snp_compound_correlation)
  
  #find_significant_cases
  ind_significant =  ((qval_pcc_snp_compound_correlation < thres_qval) & (qval_spearman_snp_compound_correlation < thres_qval))
  non_zero_count = apply(data_pca_tsne[,colnames(pcc_snp_compound_correlation)], 2, function(x){length(which(x > 0))})
  ind_val_significant = ((abs(pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres]) >= thres_corr) & (abs(spearman_snp_compound_correlation[,non_zero_count >= non_zero_count_thres]) >= thres_corr))
  ind_final_significant = ind_significant[,non_zero_count >= non_zero_count_thres] & ind_val_significant
  


  #find the significant pairs
  significant_snp_compound_pairs = which(ind_final_significant, arr.ind = T)
  significant_snp_compound_df = data.frame(SNP = rownames(pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres])[significant_snp_compound_pairs[,1]],
                                           CompoundID = colnames(pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres])[significant_snp_compound_pairs[,2]],
                                           Compound = sapply(colnames(pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres])[significant_snp_compound_pairs[,2]], map_compound_name2),
                                           PCC = pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs],
                                           Spearman = spearman_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs],
                                           Pval_PCC = pval_pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs],
                                           Pval_Spearman = pval_spearman_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs],
                                           Qval_PCC = qval_pcc_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs],
                                           Qval_Spearman = qval_spearman_snp_compound_correlation[,non_zero_count >= non_zero_count_thres][significant_snp_compound_pairs]
  )
}


#Intersect with bioactivity compounds
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))){
  stop("Missing bioactivity compound data (Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv).")
}

bioactivity_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))
bioactivity_compound = unique(setdiff(bioactivity_data$DMDID,NA))
bioactivity_compound = c(bioactivity_compound, "DMD301503", "DMD301501","DMD301505") #Beta-casein, Alpha-S1-casein, Kappa-casein

significant_snp_compound_df_bioactive = significant_snp_compound_df[which(significant_snp_compound_df$CompoundID %in% bioactivity_compound),]


#Pair/SNP/Compound counts
thres_values = seq(0.4,0.8,0.05)
snp_count = c()
compound_count = c()

snp_bioactive_count = c()
bioactive_compound_count = c()

snp_compound_pair_count = c()
snp_bioactive_compound_pair_count = c()

for (cur_thres in thres_values){
  tmp_df = significant_snp_compound_df[which(abs(significant_snp_compound_df$PCC) >= cur_thres & 
                                               abs(significant_snp_compound_df$Spearman) >= cur_thres),]
  
  tmp_df_bioactive = tmp_df[which(tmp_df$CompoundID %in% bioactivity_compound),]
  
  snp_count = c(snp_count,length(unique(tmp_df$SNP)))
  compound_count = c(compound_count,length(unique(tmp_df$CompoundID)))
  
  snp_bioactive_count = c(snp_bioactive_count,length(unique(tmp_df_bioactive$SNP)))
  bioactive_compound_count = c(bioactive_compound_count,length(unique(tmp_df_bioactive$CompoundID)))
  
  snp_compound_pair_count = c(snp_compound_pair_count, nrow(tmp_df))
  snp_bioactive_compound_pair_count = c(snp_bioactive_compound_pair_count, nrow(tmp_df_bioactive))
}

df_curve_count = data.frame("Correlation Threshold" = rep(thres_values,6),
                            "Count" = c(snp_count,
                                      compound_count,
                                      snp_bioactive_count,
                                      bioactive_compound_count,
                                      snp_compound_pair_count,
                                      snp_bioactive_compound_pair_count),
                            "Type" = factor(c(rep("#SNP in SNP-Compound Pairs", length(thres_values)),
                                       rep("#Compound in SNP-Compound Pairs", length(thres_values)),
                                       rep("#SNP in SNP-Bioactive Compound Pairs", length(thres_values)),
                                       rep("#Compound in SNP-Bioactive Compound Pairs", length(thres_values)),
                                       rep("#SNP-Compound Pairs", length(thres_values)),
                                       rep("#SNP-Bioactive Compound Pairs", length(thres_values))),
                                       levels=c("#SNP in SNP-Compound Pairs",
                                                "#Compound in SNP-Compound Pairs",
                                                "#SNP in SNP-Bioactive Compound Pairs",
                                                "#Compound in SNP-Bioactive Compound Pairs",
                                                "#SNP-Compound Pairs",
                                                "#SNP-Bioactive Compound Pairs")),
                            check.names=F
                            )

library(patchwork)
library(ggplot2)
h_pairs = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#SNP-Compound Pairs","#SNP-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#SNP-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#SNP-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  xlim(c(0.4,0.8)) + 
  ylim(c(-6000,35000)) +
  xlab("") + 
  ylab("")



h_snp = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#SNP in SNP-Compound Pairs","#SNP in SNP-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#SNP in SNP-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#SNP in SNP-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  xlim(c(0.4,0.8)) + 
  ylim(c(-2500,15000)) + 
  xlab("") + 
  ylab("")


h_compound = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#Compound in SNP-Compound Pairs","#Compound in SNP-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#Compound in SNP-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#Compound in SNP-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none") + 
  xlim(c(0.4,0.8)) + 
  xlab("") + 
  ylab("") + 
  scale_y_continuous(
    limits = c(-500,3000),
    breaks = c(0, 1000, 2000, 3000), 
    labels = paste0(" ",c(0, 1000, 2000, 3000)) 
  )

h = h_pairs / h_snp / h_compound + plot_layout(heights = c(1,1,1.2))
svg(paste0(FIGURE_PATH,"/SuppFig29a_SNP_CompoundPairs.svg"), width = 8, height = 5)
print(h)
dev.off()

#"#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462"