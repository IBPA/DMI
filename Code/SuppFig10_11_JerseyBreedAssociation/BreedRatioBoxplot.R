#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))
library(stringr)

if (!file.exists(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_adfreq_matrix.rdata")) ||
    !file.exists(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv"))){
  stop("Please run CowbreedCompositionRatioEstimation.R first to estimate the cow breed composition ratio in milk samples.")
}


sample_code_mapping = read.csv(paste0(DATA_PATH,"/SampleCodeMapping.csv"), header = T, check.names = F)
estimation_results = read.csv(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv"),row.names=1)
rownames(estimation_results) = sapply(strsplit(rownames(estimation_results),"_"), function(x){
  return(paste0(x[1:2],collapse = "_"))
})

mapped_id = sapply(rownames(estimation_results),function(x){
  idx = which(sample_code_mapping$SampleCode == x)
  if (length(idx) == 0) return(NA)
  return(sample_code_mapping$`DMD ID`[idx])
})

estimation_results = estimation_results[!is.na(mapped_id),]
rownames(estimation_results) = mapped_id[!is.na(mapped_id)]

df_jersey_ratio_region = data.frame(region = sapply(rownames(estimation_results), function(x){
  return(regions[product_data$`Region Code`[which(product_data$`DMD ID` == x)]])
}),
  jersey_ratio = estimation_results$Jersey_val
)

df_jersey_ratio_region$region = factor(df_jersey_ratio_region$region,
                                        levels = rev(regions))

pval_diff = sapply(split(df_jersey_ratio_region$jersey_ratio, df_jersey_ratio_region$region), function(x){
  other_values = df_jersey_ratio_region$jersey_ratio[!df_jersey_ratio_region$jersey_ratio %in% x]
  return(t.test(x, other_values)$p.value)
})

library(ggplot2)
library(RColorBrewer)
svg(paste0(FIGURE_PATH,"/SuppFig11a_Region_JerseyBreed.svg"), width = 4, height = 4)
p = ggplot(df_jersey_ratio_region, aes(y=region, x=jersey_ratio*100, fill=region)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  labs(title="", y="", x="Jersey Breed Ratio (%)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_fill_manual(values = rev(brewer.pal(length(regions),"Set3"))) +
  xlim(c(0,50)) + scale_y_discrete(expand = expansion(add = c(1, 1.5))) + 
  geom_jitter(width = 0, height = 0.2, alpha = 0.3, size=0.8, show.legend = F) + 
  geom_vline(xintercept = mean(df_jersey_ratio_region$jersey_ratio*100), linetype="dashed", color = "red") + 
  annotate("text", x=mean(df_jersey_ratio_region$jersey_ratio*100)+10, y=length(regions)+1, label=paste0("Mean:",round(mean(df_jersey_ratio_region$jersey_ratio)*100,2),"%"), color="red", size=4) +
  annotate("text", x=50, 
           y=(1:10)-0.15, 
           label=sapply(pval_diff,function(x){if (x < 0.05) return("*"); return("")}), size=5)
print(p)
dev.off()


#Jersey ratio vs Kappa-casein conc?
cor_jersey = apply(ordered_compound_percentage, 2, function(x){
  return(cor(as.numeric(x), as.numeric(estimation_results[names(x),"Jersey_val"]), method = "pearson"))
})
spearman_jersey = apply(ordered_compound_percentage, 2, function(x){
  return(cor(as.numeric(x), as.numeric(estimation_results[names(x),"Jersey_val"]), method = "spearman"))
})

pval_cor_jersey = apply(ordered_compound_percentage, 2, function(x){
  return(cor.test(as.numeric(x), as.numeric(estimation_results[names(x),"Jersey_val"]), method = "pearson"))
})
qval_cor_jersey = p.adjust(sapply(pval_cor_jersey, function(x){x$p.value}), method = "BH")

pval_spearman_jersey = apply(ordered_compound_percentage, 2, function(x){
  return(cor.test(as.numeric(x), as.numeric(estimation_results[names(x),"Jersey_val"]), method = "spearman"))
})
qval_spearman_jersey = p.adjust(sapply(pval_spearman_jersey, function(x){x$p.value}), method = "BH")
significant_cases = which(qval_cor_jersey < 0.1 & qval_spearman_jersey < 0.1)

#Scatter plot
df_jersey_kappa_casein = data.frame(
  Percentage = ordered_compound_percentage[, "DMD301505"]*100,
  JerseyRatio = estimation_results[rownames(ordered_compound_percentage),"Jersey_val"],
  Region = factor(regions[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)

svg(paste0(FIGURE_PATH,"/SuppFig11b_Jersey_KappaCaseinCorr.svg"), width = 5, height = 4)
p_corr2 = ggplot(df_jersey_kappa_casein, aes(y=Percentage, x=JerseyRatio, alpha=Time, shape = Brand, fill = Region)) +
  geom_point(size=3, stroke=1.0) +
  theme_classic() +
  labs(title="", x="Estimated Jersey Ratio (%)", y="Kappa-casein content (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm")) + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
  guides(fill=guide_legend(override.aes=list(shape=21)), 
         alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
         shape=guide_legend(override.aes=list(alpha=1)))
print(p_corr2)
dev.off()


#Prepare the table of compounds correlated to Jersey ratio
df_compound_jersey_correlation = data.frame(
  CompoundID = names(cor_jersey[significant_cases]),
  CompoundName = sapply(names(cor_jersey[significant_cases]), map_compound_name2),
  Content = colMeans(ordered_compound_percentage_wholemilk[,names(cor_jersey[significant_cases])])*100,
  Pearson_Correlation = cor_jersey[significant_cases],
  Pearson_Pvalue = sapply(pval_cor_jersey[significant_cases], function(x){x$p.value}),
  Pearson_Qvalue = qval_cor_jersey[significant_cases],
  Spearman_Correlation = spearman_jersey[significant_cases],
  Spearman_Pvalue = sapply(pval_spearman_jersey[significant_cases], function(x){x$p.value}),
  Spearman_Qvalue = qval_spearman_jersey[significant_cases]
)

write.csv(df_compound_jersey_correlation, paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/Compound_Correlated_to_Jersey_Breed.csv"), row.names = F)
