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

load(paste0(PROCESSED_DATA_PATH,"/SuppFig10_11_JerseyBreedAssociation/breed_composition_estimation_adfreq_matrix.rdata"))

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

colnames(ad_freq) = sapply(strsplit(colnames(ad_freq),"_"), function(x){
  return(paste0(x[1:2],collapse = "_"))
})
mapped_id_ad_freq = sapply(colnames(ad_freq),function(x){
  idx = which(sample_code_mapping$SampleCode == x)
  if (length(idx) == 0) return(NA)
  return(sample_code_mapping$`DMD ID`[idx])
})
colnames(ad_freq) = mapped_id_ad_freq
ad_freq = ad_freq[,!is.na(colnames(ad_freq))]
ad_freq_subset = ad_freq[intersect_snps,rownames(estimation_results)]
matrix_to_plot = ad_freq_subset[,product_data$`DMD ID`]

#Do KNN to fill the missing values in the matrix_to_plot
#Here, we use the pre-filled data
load(paste0(DATA_PATH,"/SNPInformation/knn_impute_result_p_10000_k_5.rdata")) 

matrix_to_plot = knn_impute_result$data

#Plot the heatmap with Jersey and Holstein composition ratios
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)

load(paste0(DATA_PATH,"/8MSNPs.rdata"))
rownames(bgvd_data_intersect) = paste0(bgvd_data_intersect$UMD3_1_1_Pos,"_",bgvd_data_intersect$Alleles)
bgvd_data_intersect = bgvd_data_intersect[intersect_snps,]

SNP_impact_level = read.csv(paste0(DATA_PATH,"/SNPInformation/SNPEffImpactLevel.csv"), header = T, check.names = F)
impact_level_dict = lapply(SNP_impact_level$Impact, function(x){x})
names(impact_level_dict) = SNP_impact_level$`Effect Seq. Ontology`
bgvd_data_intersect$impact_levels = sapply(bgvd_data_intersect$ConsequenceType, function(x){impact_level_dict[[x]]})
bgvd_data_intersect$impact_levels = stringr::str_to_title(bgvd_data_intersect$impact_levels)

snps_significant = rownames(bgvd_data_intersect)[which(bgvd_data_intersect$impact_levels %in% c("High","Moderate","Low"))]

#Annotation
annotation_all = data.frame(
  Region = product_data$`Region Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(rownames(data_pca_tsne), product_data$`DMD ID`)]
)
rownames(annotation_all) = product_data$`DMD ID`[match(rownames(data_pca_tsne), product_data$`DMD ID`)]
annotation_all$Region = factor(regions[annotation_all$Region],levels=regions)
annotation_all$Time = factor(c("T1","T2","T3")[annotation_all$Time], levels=c("T1","T2","T3"))
annotation_all$Brand = factor(c("Brand A","Brand B")[annotation_all$Brand], levels=c("Brand A","Brand B"))

annotation_wholemilk = data.frame(
  Region = product_data$`Region Code`[match(rownames(data_pca_tsne[which(!rownames(data_pca_tsne) %in% c("DMD100028","DMD100030")),]), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(rownames(data_pca_tsne[which(!rownames(data_pca_tsne) %in% c("DMD100028","DMD100030")),]), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(rownames(data_pca_tsne[which(!rownames(data_pca_tsne) %in% c("DMD100028","DMD100030")),]), product_data$`DMD ID`)]
)
rownames(annotation_wholemilk) = match(rownames(data_pca_tsne[which(!rownames(data_pca_tsne) %in% c("DMD100028","DMD100030")),]), product_data$`DMD ID`)
annotation_wholemilk$Region = factor(regions[annotation_wholemilk$Region],levels=regions)
annotation_wholemilk$Time = factor(c("T1","T2","T3")[annotation_wholemilk$Time], levels=c("T1","T2","T3"))
annotation_wholemilk$Brand = factor(c("Brand A","Brand B")[annotation_wholemilk$Brand], levels=c("Brand A","Brand B"))




#range_matrix_to_plot = apply(matrix_to_plot, 1, function(x){max(x)-min(x)})
#range_matrix_to_plot[is.na(range_matrix_to_plot)] = 0

matrix_to_plot2 = sapply(1:length(regions), function(x){
  region_samples = product_data$`DMD ID`[which(product_data$`Region Code` == x)]
  return(rowMeans(matrix_to_plot[,region_samples,drop=F]))
})

#tmp_kmeans = kmeans(matrix_to_plot, 10, algorithm="Lloyd", iter.max = 1000)
#row_order_idx = order(tmp_kmeans$cluster)

mean_estimation_results = sapply(1:length(regions), function(x){
  region_samples = product_data$`DMD ID`[which(product_data$`Region Code` == x)]
  return(colMeans(estimation_results[region_samples,,drop=F]))
})
mean_estimation_results = t(mean_estimation_results)


ha = HeatmapAnnotation(df = annotation_all[colnames(matrix_to_plot),3:1],
                       col = list(
                         Brand = setNames(c("#E41A1C","#377EB8"),c("Brand A","Brand B")),
                         Time = setNames(brewer.pal(length(unique(annotation_all$Time)),"Set2"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation_all$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left")

ha2 = HeatmapAnnotation(df = data.frame(Region=factor(regions,levels=regions)),
                        col = list(
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F)



#With the curves

ha_curve = HeatmapAnnotation("Breed\nComposition\nRatio" = anno_lines(as.matrix(estimation_results[colnames(matrix_to_plot),]),
                                               gp = gpar(col = c("#6A3D9A","#B15928"), lwd=2),
                                               add_points = T,
                                               pt_gp = gpar(col = c("#6A3D9A","#B15928"), cex=0.5, pch = 21),
                                               axis = TRUE,
                                               axis_param = list(at = c(0,0.25,0.5,0.75,1.0),
                                                                 labels = c("0","0.25","0.50","0.75","1.00")),
                                               height = unit(2, "cm")),
                       show_legend = F,
                       show_annotation_name = F)


ha_curve2 = HeatmapAnnotation("Breed\nComposition\nRatio" = anno_lines(as.matrix(mean_estimation_results),
                                                                    gp = gpar(col = c("#6A3D9A","#B15928"), lwd=2),
                                                                    add_points = T,
                                                                    pt_gp = gpar(col = c("#6A3D9A","#B15928"), cex=0.5, pch = 21),
                                                                    axis = TRUE,
                                                                    axis_param = list(at = c(0,0.25,0.5,0.75,1.0),
                                                                                      labels = c("","","","","")),
                                                                    height = unit(2, "cm")),
                             show_legend = F,
                             show_annotation_name = F)


pcc_dist_func = function(x){
  print(dim(t(x)))
  as.dist(1 - cor(t(x), method="pearson"))
}


range_ref_ad_freq_diff = apply(val_data_three_breeds_matrix_numeric[intersect_snps,], 1, function(x){abs(x[2]-x[3])})
#idx_selected = which(range_ref_ad_freq_diff >= q_select)
library(fastcluster)
idx_selected = 1 : nrow(matrix_to_plot)
tmp_kmeans_selected = kmeans(matrix_to_plot[idx_selected,], 5, algorithm="Lloyd", iter.max = 1000)
cluster_results = tmp_kmeans_selected$cluster
cluster_order = hclust(dist((tmp_kmeans_selected$centers)))$order
cluster_results = factor(cluster_results, levels=cluster_order)
range_ref_ad_freq_diff = apply(val_data_three_breeds_matrix_numeric[intersect_snps,], 1, function(x){abs(x[2]-x[3])})
row_order_idx_selected = order(cluster_results,-range_ref_ad_freq_diff,val_data_three_breeds_matrix_numeric[intersect_snps,2])

q_select = quantile(range_ref_ad_freq_diff,0.95)
#q_select = -Inf
idx_selected = which(range_ref_ad_freq_diff[row_order_idx_selected] >= q_select)

svg(paste0(FIGURE_PATH,"/SuppFig10a_SNP_Heatmap.svg"), width = 8, height = 4)
#svg("Figures\\SNP_Heatmap_NoSelected.svg", width = 8, height = 4)
h= Heatmap(matrix_to_plot[row_order_idx_selected[idx_selected],],
           row_split = cluster_results[row_order_idx_selected[idx_selected]],
           name = "Allele\nFrequency",
           top_annotation = ha,
           bottom_annotation = ha_curve,
           show_row_names = F,
           show_column_names = F,
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
           use_raster = T,
           raster_by_magick = TRUE
           #clustering_distance_columns = pcc_dist_func,
           #column_split = 3
)

h2 = Heatmap(matrix_to_plot2[row_order_idx_selected[idx_selected],],
             name = "Allele\nFrequency",
             row_split = cluster_results[row_order_idx_selected[idx_selected]],
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
             use_raster = T,
             raster_by_magick = TRUE 
             #clustering_distance_columns = pcc_dist_func,
             #column_split = 3
)
ref_ad = val_data_three_breeds_matrix_numeric[rownames(matrix_to_plot),c(2,3)]
colnames(ref_ad) = c("Holstein","Jersey")
h3 = Heatmap(ref_ad[row_order_idx_selected[idx_selected],],
             name = "Allele\nFrequency",
             row_split = cluster_results[row_order_idx_selected[idx_selected]],
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
             use_raster = T,
             raster_by_magick = TRUE 
             )

draw(h+h2+h3, merge_legend = T, ht_gap=unit(0.2,"cm"))
dev.off()



range_ref_ad_freq_diff = apply(val_data_three_breeds_matrix_numeric[intersect_snps,], 1, function(x){abs(x[2]-x[3])})
idx_selected = which(range_ref_ad_freq_diff >= q_select)
library(fastcluster)
#idx_selected = 1 : nrow(matrix_to_plot)
tmp_kmeans_selected = kmeans(matrix_to_plot[idx_selected,], 5, algorithm="Lloyd", iter.max = 1000)
cluster_results = tmp_kmeans_selected$cluster
cluster_order = hclust(dist((tmp_kmeans_selected$centers)))$order
cluster_results = factor(cluster_results, levels=cluster_order)
range_ref_ad_freq_diff = apply(val_data_three_breeds_matrix_numeric[intersect_snps,], 1, function(x){abs(x[2]-x[3])})
row_order_idx_selected = order(cluster_results,-range_ref_ad_freq_diff,val_data_three_breeds_matrix_numeric[intersect_snps,2])

q_select = quantile(range_ref_ad_freq_diff,0.95)
#q_select = -Inf
idx_selected = which(range_ref_ad_freq_diff[row_order_idx_selected] >= q_select)

svg(paste0(FIGURE_PATH,"/SuppFig10b_SNP_Heatmap_SelectedSNPs.svg"), width = 8, height = 4)
#svg("Figures\\SNP_Heatmap_NoSelected.svg", width = 8, height = 4)
h= Heatmap(matrix_to_plot[row_order_idx_selected[idx_selected],],
           row_split = cluster_results[row_order_idx_selected[idx_selected]],
           name = "Allele\nFrequency",
           top_annotation = ha,
           bottom_annotation = ha_curve,
           show_row_names = F,
           show_column_names = F,
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
           use_raster = T,
           raster_by_magick = TRUE
           #clustering_distance_columns = pcc_dist_func,
           #column_split = 3
)

h2 = Heatmap(matrix_to_plot2[row_order_idx_selected[idx_selected],],
             name = "Allele\nFrequency",
             row_split = cluster_results[row_order_idx_selected[idx_selected]],
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
             use_raster = T,
             raster_by_magick = TRUE 
             #clustering_distance_columns = pcc_dist_func,
             #column_split = 3
)
ref_ad = val_data_three_breeds_matrix_numeric[rownames(matrix_to_plot),c(2,3)]
colnames(ref_ad) = c("Holstein","Jersey")
h3 = Heatmap(ref_ad[row_order_idx_selected[idx_selected],],
             name = "Allele\nFrequency",
             row_split = cluster_results[row_order_idx_selected[idx_selected]],
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
             use_raster = T,
             raster_by_magick = TRUE 
)

draw(h+h2+h3, merge_legend = T, ht_gap=unit(0.2,"cm"))
dev.off()
