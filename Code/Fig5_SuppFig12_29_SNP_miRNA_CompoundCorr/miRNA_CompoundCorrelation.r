#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\ConcentrationMatrices.rdata"))){
  stop("Missing prepared concentration matrices. Please run the script '1_ConcentrationTablePreparation.R' first.")
}

load(paste0(PROCESSED_DATA_PATH,"\\ConcentrationMatrices.rdata"))
data_matrix_miRNA = concentration_matrix_miRNA_avg[rownames(data_pca_tsne), ]

thres_corr = 0.3
non_zero_count_thres = 10
thres_qval = 0.1

cor_miRNA_compound = cor(data_matrix_miRNA, data_pca_tsne, method = "pearson")
spearman_miRNA_compound = cor(data_matrix_miRNA, data_pca_tsne, method = "spearman")

pval_cor_miRNA_compound = apply(data_matrix_miRNA,2,function(x){
  apply(data_pca_tsne,2,function(y){
    cor.test(x,y, method = "pearson")$p.value
  })
})

pval_spearman_miRNA_compound = apply(data_matrix_miRNA,2,function(x){
  apply(data_pca_tsne,2,function(y){
    cor.test(x,y, method = "spearman")$p.value
  })
})
pval_cor_miRNA_compound = t(pval_cor_miRNA_compound)
pval_spearman_miRNA_compound = t(pval_spearman_miRNA_compound)
  

qval_cor_miRNA_compound = matrix(p.adjust(as.numeric(as.matrix(pval_cor_miRNA_compound)), method = "BH"),
                                           nrow = nrow(pval_cor_miRNA_compound),
                                           ncol = ncol(pval_cor_miRNA_compound))
rownames(qval_cor_miRNA_compound) = rownames(pval_cor_miRNA_compound)
colnames(qval_cor_miRNA_compound) = colnames(pval_cor_miRNA_compound)

qval_spearman_miRNA_compound = matrix(p.adjust(as.numeric(as.matrix(pval_spearman_miRNA_compound)), method = "BH"),
                                                nrow = nrow(pval_spearman_miRNA_compound),
                                                ncol = ncol(pval_spearman_miRNA_compound))
rownames(qval_spearman_miRNA_compound) = rownames(pval_spearman_miRNA_compound)
colnames(qval_spearman_miRNA_compound) = colnames(pval_spearman_miRNA_compound)


#find_significant_cases
ind_significant =  ((qval_cor_miRNA_compound < thres_qval) & (qval_spearman_miRNA_compound < thres_qval))
non_zero_count = apply(data_pca_tsne[,colnames(cor_miRNA_compound)], 2, function(x){length(which(x > 0))})
ind_val_significant = ((abs(cor_miRNA_compound[,non_zero_count >= non_zero_count_thres]) >= thres_corr) & (abs(spearman_miRNA_compound[,non_zero_count >= non_zero_count_thres]) >= thres_corr))
ind_final_significant = ind_significant[,non_zero_count >= non_zero_count_thres] & ind_val_significant

#find the significant pairs
significant_miRNA_compound_pairs = which(ind_final_significant, arr.ind = T)
significant_miRNA_compound_df = data.frame(miRNA = rownames(cor_miRNA_compound[,non_zero_count >= non_zero_count_thres])[significant_miRNA_compound_pairs[,1]],
                                         CompoundID = colnames(cor_miRNA_compound[,non_zero_count >= non_zero_count_thres])[significant_miRNA_compound_pairs[,2]],
                                         Compound = sapply(colnames(cor_miRNA_compound[,non_zero_count >= non_zero_count_thres])[significant_miRNA_compound_pairs[,2]], map_compound_name2),
                                         PCC = cor_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs],
                                         Spearman = spearman_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs],
                                         Pval_PCC = pval_cor_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs],
                                         Pval_Spearman = pval_spearman_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs],
                                         Qval_PCC = qval_cor_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs],
                                         Qval_Spearman = qval_spearman_miRNA_compound[,non_zero_count >= non_zero_count_thres][significant_miRNA_compound_pairs]
)





#Intersect with bioactivity compounds
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))){
  stop("Missing bioactivity compound data (Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv).")
}

bioactivity_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))
bioactivity_compound = unique(setdiff(bioactivity_data$DMDID,NA))
bioactivity_compound = c(bioactivity_compound, "DMD301503", "DMD301501","DMD301505") #Beta-casein, Alpha-S1-casein, Kappa-casein
significant_miRNA_compound_df_bioactive = significant_miRNA_compound_df[which(significant_miRNA_compound_df$CompoundID %in% bioactivity_compound),]
write.csv(significant_miRNA_compound_df_bioactive, "ProcessedData/miRNA_BioactiveCompoundCorr_SignificantPairs_10_q01.csv", row.names = F)
#Pair/miRNA/Compound counts
thres_values = seq(0.35,0.75,0.05)
miRNA_count = c()
compound_count = c()

miRNA_bioactive_count = c()
bioactive_compound_count = c()

miRNA_compound_pair_count = c()
miRNA_bioactive_compound_pair_count = c()

for (cur_thres in thres_values){
  tmp_df = significant_miRNA_compound_df[which((abs(significant_miRNA_compound_df$PCC) >= cur_thres) 
                                               & (abs(significant_miRNA_compound_df$Spearman) >= cur_thres)),]
  
  tmp_df_bioactive = tmp_df[which(tmp_df$CompoundID %in% bioactivity_compound),]
  
  miRNA_count = c(miRNA_count,length(unique(tmp_df$miRNA)))
  compound_count = c(compound_count,length(unique(tmp_df$CompoundID)))
  
  miRNA_bioactive_count = c(miRNA_bioactive_count,length(unique(tmp_df_bioactive$miRNA)))
  bioactive_compound_count = c(bioactive_compound_count,length(unique(tmp_df_bioactive$CompoundID)))
  
  miRNA_compound_pair_count = c(miRNA_compound_pair_count, nrow(tmp_df))
  miRNA_bioactive_compound_pair_count = c(miRNA_bioactive_compound_pair_count, nrow(tmp_df_bioactive))
}

df_curve_count = data.frame("Correlation Threshold" = rep(thres_values,6),
                            "Count" = c(miRNA_count,
                                        compound_count,
                                        miRNA_bioactive_count,
                                        bioactive_compound_count,
                                        miRNA_compound_pair_count,
                                        miRNA_bioactive_compound_pair_count),
                            "Type" = factor(c(rep("#miRNA in miRNA-Compound Pairs", length(thres_values)),
                                              rep("#Compound in miRNA-Compound Pairs", length(thres_values)),
                                              rep("#miRNA in miRNA-Bioactive Compound Pairs", length(thres_values)),
                                              rep("#Compound in miRNA-Bioactive Compound Pairs", length(thres_values)),
                                              rep("#miRNA-Compound Pairs", length(thres_values)),
                                              rep("#miRNA-Bioactive Compound Pairs", length(thres_values))),
                                            levels=c("#miRNA in miRNA-Compound Pairs",
                                                     "#Compound in miRNA-Compound Pairs",
                                                     "#miRNA in miRNA-Bioactive Compound Pairs",
                                                     "#Compound in miRNA-Bioactive Compound Pairs",
                                                     "#miRNA-Compound Pairs",
                                                     "#miRNA-Bioactive Compound Pairs")),
                            check.names=F
)

library(ggplot2)
library(patchwork)

h_pair = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#miRNA-Compound Pairs","#miRNA-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#miRNA-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#miRNA-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  xlim(c(0.35,0.75)) + 
  ylim(c(-800,3000)) +
  xlab("") + 
  ylab("")


h_miRNA = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#miRNA in miRNA-Compound Pairs","#miRNA in miRNA-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#miRNA in miRNA-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#miRNA in miRNA-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  xlim(c(0.35,0.75)) + 
  ylim(c(-80,300)) + 
  xlab("") + 
  ylab("")


h_compound = ggplot(df_curve_count[which(df_curve_count$Type %in% c("#Compound in miRNA-Compound Pairs","#Compound in miRNA-Bioactive Compound Pairs")),]) + 
  geom_line(aes(x=`Correlation Threshold`,y=Count,color=Type)) +
  geom_point(aes(x=`Correlation Threshold`,y=Count,fill=Type),shape=21) +
  theme_classic() + 
  scale_color_manual(values=c("#80B1D3","#8DD3C7")) +
  scale_fill_manual(values=c("#80B1D3","#8DD3C7")) +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#Compound in miRNA-Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=-1.5, size=3, color="#215874") +
  geom_text(data=df_curve_count[which(df_curve_count$Type %in% c("#Compound in miRNA-Bioactive Compound Pairs")),],
            aes(x=`Correlation Threshold`,y=Count,label=Count), vjust=2, size=3, color= "#086A5F") +
  theme(legend.position = "none") + 
  xlim(c(0.35,0.75)) + 
  xlab("") + 
  ylab("") + 
  scale_y_continuous(
    limits = c(-150,800)
    #breaks = c(0, 500, 1000, 1500, 2000, 2500), 
    #labels = paste0(" ",c(0, 500, 1000, 1500, 2000, 2500)) 
  )

h = h_pair / h_miRNA / h_compound + plot_layout(heights = c(1,1,1.2))
svg(paste0(FIGURE_PATH,"/SuppFig29b_miRNA_CompoundPairs.svg"), width = 8, height = 5)
print(h)
dev.off()






