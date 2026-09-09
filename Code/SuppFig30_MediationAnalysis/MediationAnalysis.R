#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

library(stringr)
library(lavaan)


if (!file.exists(paste0(DATA_PATH,"/SampleCodeMapping.csv"))){
  stop("Missing required sample mapping table (SampleCodeMapping.csv).")
}
if (!file.exists(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv"))){
  stop("Missing required breed composition ratio estimation results (SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_result.csv).")
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

use_prepared_mediation_results = T
if (use_prepared_mediation_results){
  if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\SuppFig30_MediationAnalysis\\lavaan_analysis_result_added_snp_metadata_manual.csv"))){
    stop("Missing required prepared mediation analysis results (lavaan_analysis_result_added_snp_metadata_manual.csv).")
  }
  
  lavaan_results = read.csv(paste0(PROCESSED_DATA_PATH,"\\SuppFig30_MediationAnalysis\\lavaan_analysis_result_added_snp_metadata_manual.csv"),row.names = 1)
}else{
  
  
  if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\SuppFig30_MediationAnalysis\\ad_matrix_snp_significant.rdata"))){
    stop("Missing required allele frequency matrix of the selected SNPs (ad_matrix_snp_significant.rdata).")
  }
  load(paste0(PROCESSED_DATA_PATH,"\\SuppFig30_MediationAnalysis\\ad_matrix_snp_significant.rdata"))
  
  if (!file.exists(paste0(PROCESSED_DATA_PATH,"/Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv"))){
    stop("Missing required SNP - Compound pairs (Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv).")
  }
  if (!file.exists(paste0(PROCESSED_DATA_PATH,"/Fig3_SuppFig6_17_Bioactivity/Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))){
    stop("Missing required target bioactivity compound table (Fig3_SuppFig6_17_Bioactivity/Plot_Efficacy_Matrix_Summary_ManualLabel.csv).")
  }
  
  snp_correlates_compound = read.csv(paste0(PROCESSED_DATA_PATH,"/Fig5_SuppFig12_29_SNP_Compound_Corr/SNP_CompoundCorr_SignificantPairs_10_q01.csv"))
  bioactivity_data = read.csv(paste0(PROCESSED_DATA_PATH,"/Fig3_SuppFig6_17_Bioactivity/Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))
  bioactivity_compound = unique(setdiff(bioactivity_data$DMDID,NA))
  bioactivity_compound = c(bioactivity_compound, "DMD301503", "DMD301501","DMD301505") #Beta-casein, Alpha-S1-casein, Kappa-casein
  
  snp_correlates_compound_selected = snp_correlates_compound[snp_correlates_compound$CompoundID %in% bioactivity_compound,]
  
  jersey_prop = estimation_results[rownames(data_pca_tsne),"Jersey_val"]
  
  lavaan_results = lapply(1:nrow(snp_correlates_compound_selected), function(i){
    print(i)
    cur_compound = data_pca_tsne[,snp_correlates_compound_selected$CompoundID[i]]
    cur_snp = ad_matrix_significant[snp_correlates_compound_selected$SNP[i], rownames(data_pca_tsne)]
    df = data.frame(
      cur_snp = cur_snp,
      jersey_prop = jersey_prop,
      y = cur_compound
    )
    
    model_txt <- '
  jersey_prop ~ a*cur_snp
  y      ~ b*jersey_prop + cprime*cur_snp
  indirect := a*b
  direct   := cprime
  total    := cprime + (a*b)
  prop_med := (a*b)/total
'
    fit <- sem(model_txt, data=df, se="bootstrap", bootstrap=2000)
    
    pe <- parameterEstimates(fit, ci = TRUE, level = 0.95, standardized = FALSE)
    
    # helper to pull a row by label
    get_by_label <- function(lbl) {
      r <- pe[pe$label == lbl, , drop = FALSE]
      if (nrow(r) == 0) return(NULL)
      r[1, ]
    }
    
    a_row <- get_by_label("a")
    b_row <- get_by_label("b")
    c_row <- get_by_label("cprime")
    ind_row <- pe[pe$lhs == "indirect", , drop = FALSE]
    dir_row <- pe[pe$lhs == "direct", , drop = FALSE]
    tot_row <- pe[pe$lhs == "total", , drop = FALSE]
    pm_row  <- pe[pe$lhs == "prop_med", , drop = FALSE]
    
    # Use CI exclusion of 0 as the robust “significant” rule for mediation
    sig_ci <- function(r) {
      if (is.null(r) || nrow(r) == 0) return(NA)
      (r$ci.lower[1] > 0) || (r$ci.upper[1] < 0)
    }
    
    out <- data.frame(
      a = a_row$est, a_ci_lo = a_row$ci.lower, a_ci_hi = a_row$ci.upper,
      b = b_row$est, b_ci_lo = b_row$ci.lower, b_ci_hi = b_row$ci.upper,
      cprime = c_row$est, cprime_ci_lo = c_row$ci.lower, cprime_ci_hi = c_row$ci.upper,
      
      indirect = ind_row$est[1], ind_ci_lo = ind_row$ci.lower[1], ind_ci_hi = ind_row$ci.upper[1],
      direct   = dir_row$est[1], dir_ci_lo = dir_row$ci.lower[1], dir_ci_hi = dir_row$ci.upper[1],
      total    = tot_row$est[1], tot_ci_lo = tot_row$ci.lower[1], tot_ci_hi = tot_row$ci.upper[1],
      prop_med = pm_row$est[1],  pm_ci_lo  = pm_row$ci.lower[1],  pm_ci_hi  = pm_row$ci.upper[1],
      
      sig_indirect = sig_ci(ind_row),
      sig_direct   = sig_ci(dir_row),
      stringsAsFactors = FALSE
    )
    
    # classify into the 4 buckets
    out$class <- with(out, ifelse(sig_indirect & !sig_direct, "Mediated-only",
                                  ifelse(!sig_indirect & sig_direct, "Direct-only",
                                         ifelse(sig_indirect & sig_direct, "Partial mediation",
                                                "Neither"))))
    list(fit = fit, out = out)
  })
  
  lavaan_results_backup = lavaan_results
  lavaan_results = do.call(rbind, lapply(lavaan_results, function(x){x$out}))
  lavaan_results = cbind(snp_correlates_compound_selected, lavaan_results)
}




jersey_correlates_compound = read.csv(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/Compound_Correlated_to_Jersey_Breed.csv"))
jersey_correlates_bioactive_compound = jersey_correlates_compound[jersey_correlates_compound$CompoundID %in% bioactivity_compound,]

lavaan_results_mediate_only = lavaan_results[which(lavaan_results$class == "Mediated-only"),]
lavaan_results_mediate_only = lavaan_results_mediate_only[which(abs(lavaan_results_mediate_only$PCC) > 0.5 & abs(lavaan_results_mediate_only$Spearman) > 0.5 & lavaan_results_mediate_only$impact_levels %in% c("MODERATE","HIGH")),]

#stacked bar for the count of the class
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidygraph)
library(ggraph)
library(colorspace)

df_count_class = data.frame(count = as.numeric(table(lavaan_results$class)),
                            class = names(table(lavaan_results$class)))
df_count_class$class = factor(df_count_class$class, levels = c("Mediated-only","Direct-only","Partial mediation","Neither"))

p_stacked = ggplot(df_count_class, aes(x = "", y = count, fill = class)) +
  geom_bar(position="stack", stat = "identity", width = 1) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "none") +
  #geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("Mediated-only" = "#E41A1C", "Direct-only" = "#377EB8", "Partial mediation" = "#7570B3", "Neither" = "grey"))

svg(paste0(FIGURE_PATH,"/SuppFig30a_StackedBar_LavaanClass.svg"), width = 1, height = 5)
print(p_stacked)
dev.off()


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

matrix_to_plot = ad_matrix_significant[lavaan_results_mediate_only$SNP, rownames(data_pca_tsne)]
rownames(matrix_to_plot) = paste0(convert_snp_name2(rownames(matrix_to_plot)),"(",lavaan_results_mediate_only$Gene,")")

matrix_to_plot2 = sapply(1:length(regions), function(x){
  region_samples = product_data$`DMD ID`[which(product_data$`Region Code` == x)]
  return(rowMeans(matrix_to_plot[,region_samples,drop=F]))
})

mean_estimation_results = sapply(1:length(regions), function(x){
  region_samples = product_data$`DMD ID`[which(product_data$`Region Code` == x)]
  return(colMeans(estimation_results[region_samples,,drop=F]))
})
mean_estimation_results = t(mean_estimation_results)

matrix_to_plot_compound = t(ordered_compound_percentage[rownames(data_pca_tsne),unique(lavaan_results_mediate_only$CompoundID)])
rownames(matrix_to_plot_compound) = sapply(rownames(matrix_to_plot_compound),map_compound_name)
matrix_to_plot_compound = t(apply(matrix_to_plot_compound, 1, function(x){
  return((x - mean(x))/mean(x) * 100)
}))

matrix_to_plot_compound2 = sapply(1:length(regions), function(x){
  region_samples = product_data$`DMD ID`[which(product_data$`Region Code` == x)]
  return(rowMeans(matrix_to_plot_compound[,region_samples,drop=F]))
})


annotation = data.frame(
  Region = product_data$`Region Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)]
)
rownames(annotation) = rownames(data_pca_tsne)
annotation$Brand = factor(paste0("Brand ",annotation$Region, c("A","B")[annotation$Brand]), 
                          levels=paste0("Brand ",sapply(1:10,function(x){rep(x,2)}), rep(c("A","B"),10)))
annotation$Region = factor(regions[annotation$Region],levels=regions)
annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))

ha = HeatmapAnnotation(df = annotation[colnames(matrix_to_plot),3:1],
                       col = list(
                         Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,desaturate(x,c(0.8)))})),
                                          levels(annotation$Brand)),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set2"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize=8),
                       simple_anno_size = unit(0.3,"cm"))

ha2 = HeatmapAnnotation(df = data.frame(Region=factor(regions,levels=regions)),
                        col = list(
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F,
                        simple_anno_size = unit(0.3,"cm"))


#With the curves
ha_curve = HeatmapAnnotation(
  "Glutamic Acid" = anno_simple(matrix_to_plot_compound["Glutamic acid",],
                                col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                                height = unit(0.3, "cm")),
  "Kappa-casein" = anno_simple(matrix_to_plot_compound["Kappa-casein",],
                               col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                               height = unit(0.3, "cm")),
  "Breed\nComposition\nRatio" = anno_lines(as.matrix(estimation_results[colnames(matrix_to_plot),]),
                                           gp = gpar(col = c("#6A3D9A","#B15928"), lwd=1),
                                           add_points = F,
                                           pt_gp = gpar(col = c("#6A3D9A","#B15928"), cex=0.5, pch = 21),
                                           axis = TRUE,
                                           axis_param = list(at = c(0,1.0),
                                                             labels = c("0","1.00")),
                                           height = unit(0.6, "cm")),
  show_legend = F,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontsize=8),
  annotation_name_side = "left",
  annotation_name_rot = 0,
  gap = unit(c(0,1), "mm"))


ha_curve2 = HeatmapAnnotation(
  "Glutamic Acid" = anno_simple(matrix_to_plot_compound2["Glutamic acid",],
                                col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                                height = unit(0.3, "cm")),
  "Kappa-casein" = anno_simple(matrix_to_plot_compound2["Kappa-casein",],
                               col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                               height = unit(0.3, "cm")),
  "Breed\nComposition\nRatio" = anno_lines(as.matrix(mean_estimation_results),
                                           gp = gpar(col = c("#6A3D9A","#B15928"), lwd=1),
                                           add_points = F,
                                           pt_gp = gpar(col = c("#6A3D9A","#B15928"), cex=0.5, pch = 21),
                                           axis = TRUE,
                                           axis_param = list(at = c(0,1.0),
                                                             labels = c("","")),
                                           height = unit(0.6, "cm")),
  show_legend = F,
  show_annotation_name = F,
  gap = unit(c(0,1), "mm"))

ref_ad = lavaan_results_mediate_only[,c("Holstein_val","Jersey_val")]
colnames(ref_ad) = c("Holstein","Jersey")
rownames(ref_ad) = rownames(matrix_to_plot)


row_order = hclust(dist(matrix_to_plot))$order

h= Heatmap(matrix_to_plot[row_order,],
           name = "Allele\nFrequency",
           top_annotation = ha,
           bottom_annotation = ha_curve,
           show_row_names = T,
           show_column_names = F,
           row_names_side="left",
           row_names_gp = gpar(fontsize=7),
           cluster_rows = F,
           #cluster_row_slices = F,
           cluster_columns = T,
           show_row_dend = F,
           show_column_dend = F,
           row_title = " ",
           #row_gap = unit(5, "mm"),
           column_title = "(All 60 Milk Samples)",
           column_title_side = "bottom",
           column_title_gp = gpar(fontsize=10),
           col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
           heatmap_legend_param = list(at=c(0,0.25,0.5,0.75,1.00),labels=c("0","0.25","0.50","0.75","1.00")),
           use_raster = F,
           raster_by_magick = F
           #clustering_distance_columns = pcc_dist_func,
           #column_split = 3
)

h2 = Heatmap(matrix_to_plot2[row_order,],
             name = "Allele\nFrequency",
             top_annotation = ha2,
             bottom_annotation = ha_curve2,
             row_names_side="left",
             show_row_names = F,
             show_column_names = F,
             cluster_rows = F,
             #cluster_row_slices = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             row_title = " ",
             #row_gap = unit(5, "mm"),
             column_title = "(Region Avg.)",
             column_title_side="bottom",
             column_title_gp = gpar(fontsize=10),
             col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0,0.25,0.5,0.75,1.00),labels=c("0","0.25","0.50","0.75","1.00")),
             use_raster = F,
             raster_by_magick = F 
             #clustering_distance_columns = pcc_dist_func,
             #column_split = 3
)

h3 = Heatmap(ref_ad[row_order,],
             name = "Allele\nFrequency",
             column_names_side="top",
             show_row_names = F,
             show_column_names = F,
             cluster_rows = F,
             cluster_columns = F,
             show_row_dend = F,
             show_column_dend = F,
             row_title = " ",
             col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0,0.25,0.5,0.75,1.00),labels=c("0","0.25","0.50","0.75","1.00")),
             use_raster = F,
             raster_by_magick = F 
)

svg(paste0(FIGURE_PATH,"/SuppFig30b_heatmap.svg"),width=6,height=3)
draw(h+h2+h3, merge_legend = T, ht_gap=unit(0.2,"cm"))
dev.off()


