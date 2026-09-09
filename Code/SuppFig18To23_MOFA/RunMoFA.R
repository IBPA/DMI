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
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\SuppFig18To23_MOFA\\ad_matrix_snp_medium.rdata"))){
  stop("Missing required allele frequency matrix of the selected SNPs (ad_matrix_snp_medium.rdata).")
}
load(paste0(PROCESSED_DATA_PATH,"/ConcentrationMatrices.rdata"))
load(paste0(PROCESSED_DATA_PATH,"\\SuppFig18To23_MOFA\\ad_matrix_snp_medium.rdata"))

data_matrix_miRNA = concentration_matrix_miRNA_avg[rownames(data_pca_tsne),]

library(data.table)
library(MOFA2)

mol_classes = sapply(colnames(ordered_compound_percentage), molecular_classification)

lipid_matrix = t(as.matrix(ordered_compound_percentage[,mol_classes=="Lipids"]))
metabolites_matrix = t(as.matrix(ordered_compound_percentage[,mol_classes=="Metabolites"]))
oligosaccharides_matrix = t(as.matrix(ordered_compound_percentage[,mol_classes=="Oligosaccharides"]))
peptides_matrix = t(as.matrix(ordered_compound_percentage[,mol_classes=="Peptides"]))
proteins_matrix = t(as.matrix(ordered_compound_percentage[,mol_classes=="Proteins"]))

#ad_matrix_snp_medium_normalized = apply(ad_matrix_snp_medium, 2, function(x){x/sum(x)})
#data_matrix_miRNA_normalized = apply(data_matrix_miRNA, 1, function(x){x/sum(x)})

MOFAobject = create_mofa(list(SNP = ad_matrix_snp_medium[,rownames(ordered_compound_percentage)],
                              miRNA = t(data_matrix_miRNA)[,rownames(ordered_compound_percentage)],
                              Proteins= proteins_matrix,
                              Lipids = lipid_matrix,
                              Metabolites = metabolites_matrix,
                              Oligosaccharides = oligosaccharides_matrix,
                              Peptides = peptides_matrix
                              ))

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
data_opts$center_groups = T
data_opts$scale_groups = F
data_opts$scale_views = F
data_opts$use_float32 = F
model_opts <- get_default_model_options(MOFAobject)
#model_opts$num_factors <- 3
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter = 10000
train_opts$convergence_mode = "slow"
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(tempdir(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)

idx_product = sapply(rownames(ordered_compound_percentage), function(x){which(product_data$`DMD ID` == x)})
metadata = data.frame(sample = product_data$`DMD ID`[idx_product],
                      Region = factor(regions[product_data$`Region Code`[idx_product]], levels=regions),
                      Brand = paste0("B",product_data$`Brand Code`[idx_product]),
                      Time = factor(paste0("T",product_data$`Time Point Code`[idx_product]),levels=c("T1","T2","T3")))
samples_metadata(MOFAobject.trained) = metadata


library(RColorBrewer)
library(ggplot2)

library(patchwork)
svg(paste0(FIGURE_PATH,"/SuppFig19_MOFASummary.svg"), width = 8, height = 10)
p0a = plot_variance_explained(MOFAobject.trained, plot_total = TRUE)[[1]] + 
  theme(legend.position = "left",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.margin = margin(t = 0, r = -30, b = 0, l = 0)) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 5)) + 
  labs(fill = "Var (%)")

p0b = plot_variance_explained(MOFAobject.trained, plot_total = TRUE)[[2]] + 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10)) + 
  ylab("Total Variance\nExplained (%)") 

scale_abs_weights_f1 = get_weights(MOFAobject.trained, views = "all", factors = 1, abs=T, scale=T)
df_boxplot_weights_views_f1 = data.frame(
  View = rep(names(scale_abs_weights_f1), sapply(scale_abs_weights_f1, nrow)),
  Abs_Weight = unlist(lapply(scale_abs_weights_f1, function(x){x[,1]}))
)
df_boxplot_weights_views_f1$View = factor(df_boxplot_weights_views_f1$View, levels=names(scale_abs_weights_f1))
p1 = ggplot(df_boxplot_weights_views_f1, aes(x = View, y = Abs_Weight, fill = View)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10)) +
  ylab("Abs. Scaled\nWeight (Factor 1)") +
  xlab("Omics Layers") + 
  scale_fill_manual(values = c("#E41A1C","#377EB8","#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462")) 

scale_abs_weights_f2 = get_weights(MOFAobject.trained, views = "all", factors = 2, abs=T, scale=T)
df_boxplot_weights_views_f2 = data.frame(
  View = rep(names(scale_abs_weights_f2), sapply(scale_abs_weights_f2, nrow)),
  Abs_Weight = unlist(lapply(scale_abs_weights_f2, function(x){x[,1]}))
)
df_boxplot_weights_views_f2$View = factor(df_boxplot_weights_views_f2$View, levels=names(scale_abs_weights_f2))
p2 = ggplot(df_boxplot_weights_views_f2, aes(x = View, y = Abs_Weight, fill = View)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10)) +
  ylab("Abs. Scaled\nWeight (Factor 2)") +
  xlab("Omics Layers") + 
  scale_fill_manual(values = c("#E41A1C","#377EB8","#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"))

p0a /p0b / p1 / p2 + plot_layout(heights = c(1.5, 1, 1, 1)) & theme(plot.margin = margin(30, 5, 0, 5))
dev.off()


#Loading of factor 1 (or 2)
weights = get_weights(MOFAobject.trained, views = "all", factors = 2)

plot_weight = function(weights){
  df_snp = data.frame(Feature = rownames(weights$SNP), Weight = weights$SNP[,1], View = "SNP")
  df_snp = df_snp[order(abs(df_snp$Weight), decreasing = TRUE),][1:20,]
  df_snp$gt0 = factor(df_snp$Weight > 0, levels=c(T,F))
  df_snp$Feature = gsub("GK0000","",df_snp$Feature)
  df_snp$Feature = gsub("GK00000","",df_snp$Feature)
  df_snp$Feature = gsub(".2","",df_snp$Feature,fixed=T)
  df_snp$Feature = sapply(df_snp$Feature, function(x){
    tmp = strsplit(x, "_")[[1]]
    return(paste0(tmp[1], ":", tmp[2], ":", tmp[3], "/", tmp[4]))
  })
  df_snp$Feature = factor(df_snp$Feature, levels=rev(df_snp$Feature))
  p_snp = ggplot(df_snp, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 SNPs") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "grey")
  
  df_miRNA = data.frame(Feature = rownames(weights$miRNA), Weight = weights$miRNA[,1], View = "miRNA")
  df_miRNA = df_miRNA[order(abs(df_miRNA$Weight), decreasing = TRUE),][1:20,]
  df_miRNA$gt0 = factor(df_miRNA$Weight > 0, levels=c(T,F))
  df_miRNA$Feature = sapply(df_miRNA$Feature, function(x){idx=which(mol_data$`DMD ID`==x); mol_data$`External Database IDs`[idx]})
  df_miRNA$Feature = gsub('[]{}[\"]', "",df_miRNA$Feature)
  df_miRNA$Feature = factor(df_miRNA$Feature, levels=rev(df_miRNA$Feature))
  
  p_miRNA = ggplot(df_miRNA, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 miRNAs") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  df_protein = data.frame(Feature = rownames(weights$Proteins), Weight = weights$Proteins[,1], View = "Proteins")
  df_protein = df_protein[order(abs(df_protein$Weight), decreasing = TRUE),][1:20,]
  df_protein$gt0 = factor(df_protein$Weight > 0, levels=c(T,F))
  df_protein$Feature = sapply(df_protein$Feature, map_compound_name2)
  df_protein$Feature = gsub("\n","",df_protein$Feature)
  df_protein$Feature = factor(df_protein$Feature, levels=rev(df_protein$Feature))
  
  p_protein = ggplot(df_protein, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 Proteins") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  df_lipid = data.frame(Feature = rownames(weights$Lipids), Weight = weights$Lipids[,1], View = "Lipids")
  df_lipid = df_lipid[order(abs(df_lipid$Weight), decreasing = TRUE),][1:20,]
  df_lipid = df_lipid[df_lipid$Weight != 0,]
  df_lipid$gt0 = factor(df_lipid$Weight > 0, levels=c(T,F))
  df_lipid$Feature = sapply(df_lipid$Feature, map_compound_name)
  df_lipid$Feature = factor(df_lipid$Feature, levels=rev(df_lipid$Feature))
  
  p_lipid = ggplot(df_lipid, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 Lipids") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  df_metabolites = data.frame(Feature = rownames(weights$Metabolites), Weight = weights$Metabolites[,1], View = "Metabolites")
  df_metabolites = df_metabolites[order(abs(df_metabolites$Weight), decreasing = TRUE),][1:20,]
  df_metabolites = df_metabolites[df_metabolites$Weight != 0,]
  df_metabolites$gt0 = factor(df_metabolites$Weight > 0, levels=c(T,F))
  df_metabolites$Feature = sapply(df_metabolites$Feature, map_compound_name)
  df_metabolites$Feature = factor(df_metabolites$Feature, levels=rev(df_metabolites$Feature))
  
  p_metabolites = ggplot(df_metabolites, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 Metabolites") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  df_oligosaccharides = data.frame(Feature = rownames(weights$Oligosaccharides), Weight = weights$Oligosaccharides[,1], View = "Oligosaccharides")
  df_oligosaccharides = df_oligosaccharides[order(abs(df_oligosaccharides$Weight), decreasing = TRUE),][1:20,]
  df_oligosaccharides = df_oligosaccharides[df_oligosaccharides$Weight != 0,]
  df_oligosaccharides$gt0 = factor(df_oligosaccharides$Weight > 0, levels=c(T,F))
  df_oligosaccharides$Feature = sapply(df_oligosaccharides$Feature, map_compound_name)
  df_oligosaccharides$Feature = gsub("' -","'-",df_oligosaccharides$Feature, fixed=T)
  df_oligosaccharides$Feature = factor(df_oligosaccharides$Feature, levels=rev(df_oligosaccharides$Feature))
  
  p_oligosaccharides = ggplot(df_oligosaccharides, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size=8)) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 Oligosaccharides") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  
  df_peptides = data.frame(Feature = rownames(weights$Peptides), Weight = weights$Peptides[,1], View = "Peptides")
  df_peptides = df_peptides[order(abs(df_peptides$Weight), decreasing = TRUE),][1:20,]
  df_peptides = df_peptides[df_peptides$Weight != 0,]
  df_peptides$gt0 = factor(df_peptides$Weight > 0, levels=c(T,F))
  df_peptides$Feature = sapply(df_peptides$Feature, map_compound_name2)
  df_peptides$Feature = gsub("Peptide:","",df_peptides$Feature,fixed=T)
  df_peptides$Feature = factor(df_peptides$Feature, levels=rev(df_peptides$Feature))
  
  p_peptides = ggplot(df_peptides, aes(y = Feature, x = Weight, fill = gt0)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
    ylab("Top 20 Peptides") +
    xlab("Weight") + 
    geom_vline(xintercept = 0, color = "black")
  
  
  p = (p_snp | p_miRNA) / 
    (p_metabolites | p_oligosaccharides) / 
    free(p_protein) +
    free(p_lipid) +
    free(p_peptides) + 
    plot_layout(heights = c(1, 1, 1, 1, 1)) & 
    theme(plot.margin = margin(5, 5, 5, 5))
}
svg(paste0(FIGURE_PATH,"/SuppFig21_Weights_Factor1.svg"), width = 9, height = 12)
weights = get_weights(MOFAobject.trained, views = "all", factors = 1)
p = plot_weight(weights)
print(p)
dev.off()

svg(paste0(FIGURE_PATH,"/SuppFig23_Weights_Factor2.svg"), width = 9, height = 12)
weights = get_weights(MOFAobject.trained, views = "all", factors = 2)
p = plot_weight(weights)
print(p)
dev.off()


#Scatter plot for factor 1 and 2
df <- plot_factors(MOFAobject.trained, factors = 1:2, color_by="Region", shape_by="Brand", return_data = TRUE)
df$Time = metadata$Time
p1 = ggplot(df, aes(x = x, y = y, shape=shape_by, alpha=Time)) +
  geom_point(size = 3,aes(fill=color_by)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) +
  scale_fill_manual(values = brewer.pal(10,"Set3")) + 
  scale_shape_manual(values = c(21, 24)) +
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("Factor 1") +
  ylab("Factor 2")
print(p1)

df <- plot_factors(MOFAobject.trained, factors = 2:3, color_by="Region", shape_by="Brand", return_data = TRUE)
df$Time = metadata$Time
p2 = ggplot(df, aes(x = x, y = y, fill=color_by, shape=shape_by, alpha=Time)) +
  geom_point(size = 3) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) +
  scale_fill_manual(values = brewer.pal(10,"Set3")) + 
  scale_shape_manual(values = c(21, 24)) +
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("Factor 2") +
  ylab("Factor 3")
print(p2)

p = p1 + p2
svg(paste0(FIGURE_PATH,"/SuppFig18_factor_values.svg"), width = 6, height = 3)
print(p)
dev.off()




MOFAobject = create_mofa(list(SNP = ad_matrix_snp_medium[,rownames(ordered_compound_percentage)],
                              miRNA = t(data_matrix_miRNA)[,rownames(ordered_compound_percentage)],
                              Proteins= proteins_matrix,
                              Lipids = lipid_matrix,
                              Metabolites = metabolites_matrix,
                              Oligosaccharides = oligosaccharides_matrix,
                              Peptides = peptides_matrix
))

library(ComplexHeatmap)
library(colorspace)
library(circlize)
#Actual heatmap plot
#Features in Factor N

plot_weight_heatmap = function(weights){
  weights_snp = weights$SNP[order(abs(weights$SNP[,1]),decreasing=T)[1:10],1]
  matrix_snp = ad_matrix_snp_medium[names(weights_snp),rownames(ordered_compound_percentage)]
  matrix_snp = t(apply(matrix_snp, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  
  bar_colors_snp = ifelse(weights_snp >= 0, "#377EB8", "#E41A1C")
  bar_snp = rowAnnotation(
    "weights" = anno_barplot(weights_snp, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_snp),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90
  )
  
  weights_miRNA = weights$miRNA[order(abs(weights$miRNA[,1]),decreasing=T)[1:10],1]
  matrix_miRNA = t(data_matrix_miRNA)[names(weights_miRNA),rownames(ordered_compound_percentage)]
  matrix_miRNA = t(apply(matrix_miRNA, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_miRNA) = sapply(rownames(matrix_miRNA), function(x){idx=which(mol_data$`DMD ID`==x); mol_data$`External Database IDs`[idx]})
  rownames(matrix_miRNA) = gsub('[]{}[\"]', "",rownames(matrix_miRNA))
  
  bar_colors_miRNA = ifelse(weights_miRNA >= 0, "#377EB8", "#E41A1C")
  bar_miRNA = rowAnnotation(
    "weights" = anno_barplot(weights_miRNA, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_miRNA),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  weights_proteins = weights$Proteins[order(abs(weights$Proteins[,1]),decreasing=T)[1:10],1]
  matrix_proteins = proteins_matrix[names(weights_proteins),rownames(ordered_compound_percentage)]
  matrix_proteins = t(apply(matrix_proteins, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_proteins) = gsub("\n"," ",sapply(rownames(matrix_proteins),map_compound_name2))
  
  bar_colors_proteins = ifelse(weights_proteins >= 0, "#377EB8", "#E41A1C")
  bar_proteins = rowAnnotation(
    "weights" = anno_barplot(weights_proteins, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_proteins),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  weights_lipids = weights$Lipids[order(abs(weights$Lipids[,1]),decreasing=T)[1:10],1]
  matrix_lipids = lipid_matrix[names(weights_lipids),rownames(ordered_compound_percentage)]
  matrix_lipids = t(apply(matrix_lipids, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_lipids) = sapply(rownames(matrix_lipids),map_compound_name)
  
  bar_colors_lipids = ifelse(weights_lipids >= 0, "#377EB8", "#E41A1C")
  bar_lipids = rowAnnotation(
    "weights" = anno_barplot(weights_lipids, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_lipids),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  weights_metabolites = weights$Metabolites[order(abs(weights$Metabolites[,1]),decreasing=T)[1:10],1]
  matrix_metabolites = metabolites_matrix[names(weights_metabolites),rownames(ordered_compound_percentage)]
  matrix_metabolites = t(apply(matrix_metabolites, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_metabolites) = sapply(rownames(matrix_metabolites),map_compound_name)
  
  bar_colors_metabolites = ifelse(weights_metabolites >= 0, "#377EB8", "#E41A1C")
  bar_metabolites = rowAnnotation(
    "weights" = anno_barplot(weights_metabolites, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_metabolites),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  weights_peptides = weights$Peptides[order(abs(weights$Peptides[,1]),decreasing=T)[1:10],1]
  matrix_peptides = peptides_matrix[names(weights_peptides),rownames(ordered_compound_percentage)]
  matrix_peptides = t(apply(matrix_peptides, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_peptides) = sapply(rownames(matrix_peptides),map_compound_name2)
  
  bar_colors_peptides = ifelse(weights_peptides >= 0, "#377EB8", "#E41A1C")
  bar_peptides = rowAnnotation(
    "weights" = anno_barplot(weights_peptides, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_peptides),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  weights_oligosaccharides = weights$Oligosaccharides[order(abs(weights$Oligosaccharides[,1]),decreasing=T)[1:10],1]
  matrix_oligosaccharides = oligosaccharides_matrix[names(weights_oligosaccharides),rownames(ordered_compound_percentage)]
  matrix_oligosaccharides = t(apply(matrix_oligosaccharides, 1, function(x){
    (x - mean(x))/mean(x) * 100
  }))
  rownames(matrix_oligosaccharides) = gsub("' -","'-",sapply(rownames(matrix_oligosaccharides),map_compound_name), fixed=T)
  
  bar_colors_oligosaccharides = ifelse(weights_oligosaccharides >= 0, "#377EB8", "#E41A1C")
  bar_oligosaccharides = rowAnnotation(
    "weights" = anno_barplot(weights_oligosaccharides, 
                             width = unit(1, "cm"),
                             gp = gpar(fill = bar_colors_oligosaccharides),
                             axis_param = list(side = "top", 
                                               labels_rot = 45,
                                               gp=gpar(fontsize=7))),
    annotation_name_side="top",
    annotation_name_rot = 90,
    show_annotation_name = F
  )
  
  column_order = hclust(dist(t(rbind(matrix_snp,
                                     matrix_miRNA,
                                     matrix_proteins,
                                     matrix_lipids,
                                     matrix_metabolites,
                                     matrix_peptides,
                                     matrix_oligosaccharides))))$order
  
  annotation = data.frame(
    Region = product_data$`Region Code`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)],
    Time = product_data$`Time Point Code`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)],
    Brand = product_data$`Brand Code`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)]
  )
  rownames(annotation) = rownames(annotation)
  annotation$Brand = factor(paste0("Brand ",annotation$Region, c("A","B")[annotation$Brand]), 
                            levels=paste0("Brand ",sapply(1:10,function(x){rep(x,2)}), rep(c("A","B"),10)))
  annotation$Region = factor(regions[annotation$Region],levels=regions)
  annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))
  rownames(annotation) = rownames(ordered_compound_percentage)
  annotation = annotation[column_order,]
  ha = HeatmapAnnotation(Brand = annotation$Brand,
                         Time = annotation$Time,
                         Region = annotation$Region,
                         col = list(
                           Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,desaturate(x,c(0.8)))})),
                                            levels(annotation$Brand)),
                           Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set1"),c("T1","T2","T3")),
                           Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                         ),
                         show_legend = F,
                         annotation_name_side = "left",
                         annotation_name_rot = c(0,0,0),
                         annotation_name_gp = gpar(fontsize=9))
  
  h_snp = Heatmap(matrix_snp[,column_order],
                  name = "Above/Below Average (%)",
                  top_annotation = ha,
                  right_annotation = bar_snp,
                  row_names_side="left", show_row_names = T, show_column_names = F, 
                  cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                  row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                  column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                  col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                  heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                              direction = "horizontal", title_position = "leftcenter")
  )
  
  
  h_miRNA = Heatmap(matrix_miRNA[,column_order],
                    name = "Above/Below Average (%)",
                    right_annotation = bar_miRNA, 
                    row_names_side="left", show_row_names = T, show_column_names = F, 
                    cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                    row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                    column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                    col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                    heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                direction = "horizontal", title_position = "leftcenter")
  )
  
  h_protein = Heatmap(matrix_proteins[,column_order],
                      name = "Above/Below Average (%)",
                      right_annotation = bar_proteins, 
                      row_names_side="left", show_row_names = T, show_column_names = F, 
                      cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                      row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                      column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                      col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                      heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                  direction = "horizontal", title_position = "leftcenter")
  )
  
  h_lipid = Heatmap(matrix_lipids[,column_order],
                    name = "Above/Below Average (%)",
                    right_annotation = bar_lipids, 
                    row_names_side="left", show_row_names = T, show_column_names = F, 
                    cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                    row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                    column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                    col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                    heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                direction = "horizontal", title_position = "leftcenter")
  )
  
  h_metabolite = Heatmap(matrix_metabolites[,column_order],
                         name = "Above/Below Average (%)",
                         right_annotation = bar_metabolites, 
                         row_names_side="left", show_row_names = T, show_column_names = F, 
                         cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                         row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                         column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                         col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                         heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                     direction = "horizontal", title_position = "leftcenter")
  )
  
  h_peptide = Heatmap(matrix_peptides[,column_order],
                      name = "Above/Below Average (%)",
                      right_annotation = bar_peptides, 
                      row_names_side="left", show_row_names = T, show_column_names = F, 
                      cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                      row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                      column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                      col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                      heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                  direction = "horizontal", title_position = "leftcenter")
  )
  
  
  h_oligosaccharide = Heatmap(matrix_oligosaccharides[,column_order],
                              name = "Above/Below Average (%)",
                              right_annotation = bar_oligosaccharides, 
                              row_names_side="left", show_row_names = T, show_column_names = F, 
                              cluster_rows = F, cluster_row_slices = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                              row_title = " ",row_gap = unit(1.5, "mm"),row_names_gp = gpar(fontsize=7), row_names_max_width = unit(10, "cm"),
                              column_title = "", column_title_side = "bottom", column_title_gp = gpar(fontsize=10),
                              col = colorRamp2(c(-100, -50, 0, 50, 100), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
                              heatmap_legend_param = list(at=c(-100, -50, 0, 50, 100),labels=c("-100%","-50%","0%","50%","100%"), 
                                                          direction = "horizontal", title_position = "leftcenter")
  )
  
  
  ht_list = h_snp %v% h_miRNA %v% h_protein %v% h_lipid %v% h_metabolite %v% h_peptide %v% h_oligosaccharide
  return(ht_list)
}


svg(paste0(FIGURE_PATH,"/SuppFig20_WeightHeatmap_Factor1.svg"),width=9,height=12)
weights = get_weights(MOFAobject.trained, views = "all", factors = 1)
ht_list = plot_weight_heatmap(weights)
draw(ht_list, ht_gap = unit(8, "mm"), merge_legend=T, heatmap_legend_side="bottom", padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()

svg(paste0(FIGURE_PATH,"/SuppFig22_WeightHeatmap_Factor2.svg"),width=9,height=12)
weights = get_weights(MOFAobject.trained, views = "all", factors = 2)
ht_list = plot_weight_heatmap(weights)
draw(ht_list, ht_gap = unit(8, "mm"), merge_legend=T, heatmap_legend_side="bottom", padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()







