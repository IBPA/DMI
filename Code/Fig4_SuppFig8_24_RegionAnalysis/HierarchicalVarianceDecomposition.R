#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))){
  stop("Please run the script 'MergeSmallMoleculeAndPeptideBioactivity.r' with necessary manual labeling first to prepare the bioactivity efficacy data for plotting.")
}

library(lme4)
library(emmeans)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(colorspace)
library(maps)
library(mapdata)
library(scales)


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

ordered_compound_percentage_nonzero = ordered_compound_percentage[,colSums(ordered_compound_percentage) > 0]
result = lapply(colnames(ordered_compound_percentage_nonzero), function(x){
  df = data.frame(
    Value = ordered_compound_percentage_nonzero[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand
  )
  
  model = lme4::lmer(Value ~ 1 + (1|Region) + (1|Time) + (1|Brand:Region), data = df)
  var_comp = as.data.frame(VarCorr(model))
  var_comp$Proportion = var_comp$vcov/sum(var_comp$vcov)
  
  return(var_comp[,c("grp","vcov","sdcor","Proportion")])
})

tmp_matrix = t(sapply(result, function(x){
  x$Proportion
}))
colnames(tmp_matrix) = result[[1]]$grp
rownames(tmp_matrix) = colnames(ordered_compound_percentage_nonzero)

dominate_idx = (apply(tmp_matrix,1,which.max))
region_dominate_idx = which(dominate_idx == 2)
idx_region_ge_35 = which(tmp_matrix[, "Region"] >= 0.35)

#Heatmap of the propotion for the whole compounds
ha = HeatmapAnnotation(
  "Factor Dominance" = factor(c("Brand-Region","Region","Time","Residual")[dominate_idx], levels=c("Brand-Region","Region","Time","Residual")),
  "Region Proportion" = factor(c("No","Yes")[as.numeric(1:nrow(tmp_matrix) %in% idx_region_ge_35)+1], levels=c("Yes","No")),
  "Molecular Class" = factor(sapply(colnames(ordered_compound_percentage_nonzero), molecular_classification),levels=c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides")),
                   col = list(
                     "Factor Dominance" = setNames(c(brewer.pal(3,"Set1"),"#CCCCCC"), c("Brand-Region","Region","Time","Residual")),
                      "Region Proportion" = setNames(c("yellow","#EEEEEE"), c("Yes","No")),
                     "Molecular Class" = setNames(c("#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"), c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides"))
                      ),
                   annotation_legend_param = list("Factor Dominance" = list(nrow = 2),"Molecular Class" = list(nrow = 2)),
                   annotation_name_gp = gpar(fontsize = 10),
                   annotation_name_side = "left"
)
svg(paste0(FIGURE_PATH,"\\SuppFig24a_HierarchicalVarianceDecomposition.svg"), width=8, height=3)
h = Heatmap(t(tmp_matrix), name="Proportion of Variance",
        top_annotation = ha,
        column_split = c("Brand-Region","Region","Time","Residual")[dominate_idx],
        show_row_names = TRUE, show_column_names = F,
        row_names_gp = gpar(fontsize=10),
        cluster_rows = FALSE, cluster_columns = FALSE,
        use_raster = T,
        raster_by_magick = TRUE,
        heatmap_legend_param = list(direction = "horizontal"),
        col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")))
draw(h, merge_legend = T, heatmap_legend_side = "bottom")
dev.off()
  

sum_top_N_compound_percentage_per_dominate_class = lapply(c(10,20,50,100,nrow(tmp_matrix)), function(N){
  result = sapply(1:ncol(tmp_matrix),function(i){length(which(dominate_idx[1:N] == i))})/N
  names(result) = colnames(tmp_matrix)
  return(result)
})


df_top_compounds_sum_per_dominate_class_for_barplot = data.frame(
  percentage = 100*do.call(c,lapply(sum_top_N_compound_percentage_per_dominate_class,function(x){as.vector(x)})),
  Compound = factor(c(rep("Sum of top-10\nCompounds", times = length(sum_top_N_compound_percentage_per_dominate_class[[1]])),
                      rep("Sum of top-20\nCompounds", times = length(sum_top_N_compound_percentage_per_dominate_class[[1]])),
                      rep("Sum of top-50\nCompounds", times = length(sum_top_N_compound_percentage_per_dominate_class[[1]])),
                      rep("Sum of top-100\nCompounds", times = length(sum_top_N_compound_percentage_per_dominate_class[[1]])),
                      rep("Sum of all 5,172\nQuantified\nCompounds", times = length(sum_top_N_compound_percentage_per_dominate_class[[1]]))),
                    levels = c("Sum of all 5,172\nQuantified\nCompounds", 
                               "Sum of top-100\nCompounds", 
                               "Sum of top-50\nCompounds",
                               "Sum of top-20\nCompounds",
                               "Sum of top-10\nCompounds")),
  DominantClass = factor(rep(names(sum_top_N_compound_percentage_per_dominate_class[[1]]),  times = length(sum_top_N_compound_percentage_per_dominate_class)),
                          levels = rev(colnames(tmp_matrix)))
)

svg(paste0(FIGURE_PATH,"\\SuppFig24c_DominanceFactorRatio.svg"), width = 5, height = 4)

brand_region_label = paste0("",round(sapply(1:5, function(i){sum_top_N_compound_percentage_per_dominate_class[[i]]["Brand:Region"]*100}),1),"%")
brand_region_pos = sapply(1:5, function(i){
    sum_top_N_compound_percentage_per_dominate_class[[i]]["Brand:Region"]*100/2
})

region_label = paste0("",round(sapply(1:5, function(i){sum_top_N_compound_percentage_per_dominate_class[[i]]["Region"]*100}),1),"%")
region_pos = sapply(1:5, function(i){
  sum(sum_top_N_compound_percentage_per_dominate_class[[i]][c("Brand:Region")])*100 +
    sum_top_N_compound_percentage_per_dominate_class[[i]]["Region"]*100/2
})

residual_label = paste0("",round(sapply(1:5, function(i){sum_top_N_compound_percentage_per_dominate_class[[i]]["Residual"]*100}),1),"%")
residual_pos = sapply(1:5, function(i){
  sum(sum_top_N_compound_percentage_per_dominate_class[[i]][c("Brand:Region","Region")])*100 +
    sum_top_N_compound_percentage_per_dominate_class[[i]]["Residual"]*100/2
})

df_text = data.frame(
  x = c(brand_region_pos, region_pos, residual_pos),
  y = rep(c(5,4,3,2,1),3),
  label = c(brand_region_label,
            region_label,
            residual_label)
)
p = ggplot(df_top_compounds_sum_per_dominate_class_for_barplot, aes(y=Compound, x=percentage)) +
  geom_bar(stat="summary", fun = "mean", position = "stack", aes(fill=DominantClass), width=0.5) +
  #stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.3, position = position_stack(vjust = 0.5)) + # For error bars (mean +/- standard error)
  #geom_jitter(width = 0, height = 0.4, alpha = 0.1, size=0.2, show.legend = F) +
  theme_classic() +
  labs(title="", y="", x="Proportion of N compounds (%)", fill="") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "bottom", 
        legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=8))+
  scale_fill_manual(values = rev(c(brewer.pal(3,"Set1"),"#CCCCCC")),guide = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,20)) + 
  geom_text(data=df_text, aes(x=x, y=y, label=label), size=2.8)
print(p)
dev.off()

sum_top_N_compound_percentage_significant = lapply(c(10,20,50,100,nrow(tmp_matrix)), function(N){
  region_dominant = length(which(dominate_idx[1:N] == 2))/N
  ge_35 = length(setdiff(which(tmp_matrix[1:N, "Region"] >= 0.35),which(dominate_idx[1:N] == 2)))/N
  
  result = c(region_dominant, ge_35)
  names(result) = c("Region Dominant", "Region >= 35% but not dominant")
  return(result)
})


df_top_compounds_sum_significant_for_barplot = data.frame(
  percentage = 100*do.call(c,lapply(sum_top_N_compound_percentage_significant,function(x){as.vector(x)})),
  Compound = factor(c(rep("Sum of top-10\nCompounds", times = length(sum_top_N_compound_percentage_significant[[1]])),
                      rep("Sum of top-20\nCompounds", times = length(sum_top_N_compound_percentage_significant[[1]])),
                      rep("Sum of top-50\nCompounds", times = length(sum_top_N_compound_percentage_significant[[1]])),
                      rep("Sum of top-100\nCompounds", times = length(sum_top_N_compound_percentage_significant[[1]])),
                      rep("Sum of all 5,172\nQuantified\nCompounds", times = length(sum_top_N_compound_percentage_significant[[1]]))),
                    levels = c("Sum of all 5,172\nQuantified\nCompounds", 
                               "Sum of top-100\nCompounds", 
                               "Sum of top-50\nCompounds",
                               "Sum of top-20\nCompounds",
                               "Sum of top-10\nCompounds")),
  DominantClass = factor(rep(names(sum_top_N_compound_percentage_significant[[1]]),  times = length(sum_top_N_compound_percentage_significant)),
                         levels = rev(names(sum_top_N_compound_percentage_significant[[1]])))
)

svg(paste0(FIGURE_PATH,"\\SuppFig24d_RegionDominantRatio.svg"), width = 5, height = 4)

region_label = paste0("",round(sapply(1:5, function(i){sum_top_N_compound_percentage_significant[[i]]["Region Dominant"]*100}),1),"%")
region_pos = sapply(1:5, function(i){
  sum_top_N_compound_percentage_significant[[i]]["Region Dominant"]*100/2
})

ge35_label = paste0("",round(sapply(1:5, function(i){sum_top_N_compound_percentage_significant[[i]]["Region >= 35% but not dominant"]*100}),1),"%")
ge35_pos = sapply(1:5, function(i){
  sum(sum_top_N_compound_percentage_significant[[i]][c("Region Dominant")])*100 +
    sum_top_N_compound_percentage_significant[[i]]["Region >= 35% but not dominant"]*100/2
})

df_text = data.frame(
  x = c(region_pos, ge35_pos),
  y = rep(c(5,4,3,2,1),2),
  label = c(region_label,
            ge35_label)
)
p = ggplot(df_top_compounds_sum_significant_for_barplot, aes(y=Compound, x=percentage)) +
  geom_bar(stat="summary", fun = "mean", position = "stack", aes(fill=DominantClass), width=0.5) +
  theme_classic() +
  labs(title="", y="", x="Proportion of N compounds (%)", fill="") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "bottom", 
        legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=8))+
  scale_fill_manual(values = rev(c("#377EB8","#87BEFA")),guide = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,20)) + 
  geom_text(data=df_text, aes(x=x, y=y, label=label), size=2.8)
print(p)
dev.off()



#Additional analysis -- total protein, total lipid, total metabolite, total oligosaccharide, total peptide, total saccharide, cholestrol (DMD306416)
#Additional analysis -- top 20 compounds
#Additional analysis -- selected bioactive compounds
#DMD302040
#DMD302652
#DMD302690
#DMD302653
#DMD302709
#DMD302654
#DMD302287
#DMD302277
#DMD302257
#DMD302599
#DMD302655
#DMD302147

#VYPFPGPIP -> Beta-casein -> DMD301503
#IPI -> Kappa-casein -> DMD301505
#MAP -> Beta-casein -> DMD301503
#YFYPEL -> Alpha-S1-casein -> DMD301501
#EL -> Alpha-S1-casein -> DMD301501
#VYPFPGPIPN -> Beta-casein -> DMD301503

molecular_class = sapply(colnames(ordered_compound_percentage), molecular_classification)
sum_molecular_class = sapply(unique(molecular_class), function(mol_class){
  mol_ids = colnames(ordered_compound_percentage)[which(molecular_class[colnames(ordered_compound_percentage)] == mol_class)]
  return(rowSums(ordered_compound_percentage[, mol_ids, drop=F]))
})
sum_molecular_class = cbind(sum_molecular_class, rowSums(ordered_compound_percentage[, sapply(colnames(ordered_compound_percentage), function(x){
  idx = which(x == mol_data$`DMD ID`)
  if (grepl("[Ss]accharide",mol_data$`Molecule Classification`[idx])) return(T)
  return(F)
}), drop=F]))
colnames(sum_molecular_class)[ncol(sum_molecular_class)] = "Saccharides"

result_molecular_class = lapply(colnames(sum_molecular_class), function(x){
  df = data.frame(
    Value = sum_molecular_class[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand
  )
  
  model = lme4::lmer(Value ~ 1 + (1|Region) + (1|Time) + (1|Brand:Region), data = df)
  var_comp = as.data.frame(VarCorr(model))
  var_comp$Proportion = var_comp$vcov/sum(var_comp$vcov)
  
  return(var_comp[,c("grp","vcov","sdcor","Proportion")])
})
get_result_all = function(){
  df = data.frame(
    Value = rowSums(ordered_compound_percentage),
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand
  )
  model = lme4::lmer(Value ~ 1 + (1|Region) + (1|Time) + (1|Brand:Region), data = df)
  var_comp = as.data.frame(VarCorr(model))
  var_comp$Proportion = var_comp$vcov/sum(var_comp$vcov)
  
  return(var_comp[,c("grp","vcov","sdcor","Proportion")])
}
result_all = get_result_all()[,"Proportion"]

tmp_matrix_molecular_class = sapply(result_molecular_class, function(x){
  x$Proportion
})
tmp_matrix_molecular_class = t(tmp_matrix_molecular_class)
colnames(tmp_matrix_molecular_class) = result[[1]]$grp
rownames(tmp_matrix_molecular_class) = colnames(sum_molecular_class)


bioactive_compound_small_molecule = c("DMD302040","DMD302652","DMD302690","DMD302653","DMD302709","DMD302654","DMD302287","DMD302277","DMD302257","DMD302599","DMD302269","DMD302147")
bioactive_compound_small_molecule = c(bioactive_compound_small_molecule, c("DMD301501","DMD301503","DMD301505"))

top20_bioactive_compounds = intersect(bioactive_compound_small_molecule, colnames(ordered_compound_percentage)[1:20])
bioactive_compounds_only = setdiff(bioactive_compound_small_molecule, top20_bioactive_compounds)
top20_compound_only = setdiff(colnames(ordered_compound_percentage)[1:20], bioactive_compound_small_molecule)

deliver_plot_heatmap = rbind(tmp_matrix_molecular_class, 
                             tmp_matrix["DMD306416",,drop=F],
                             tmp_matrix[top20_bioactive_compounds,],
                             tmp_matrix[top20_compound_only,],
                             tmp_matrix[bioactive_compounds_only,])
deliver_column_split = c(rep("Molecular Class\n and Nutrition Facts",nrow(tmp_matrix_molecular_class)+1), 
                         rep("Bioactive Compounds\n(in top 20 compounds)",length(top20_bioactive_compounds)), 
                         rep("Other\nTop 20 Compounds",length(top20_compound_only)), 
                         rep("Other\nBioactive Compounds",length(bioactive_compounds_only)))
deliver_column_split = factor(deliver_column_split, levels=c("Molecular Class\n and Nutrition Facts", "Bioactive Compounds\n(in top 20 compounds)", "Other\nTop 20 Compounds", "Other\nBioactive Compounds"))

#Heatmap of the propotion for the whole compounds
ha = HeatmapAnnotation(
  "Factor Dominance" = factor(c("Brand-Region","Region","Time","Residual")[apply(deliver_plot_heatmap,1,which.max)], levels=c("Brand-Region","Region","Time","Residual")),
  "Region Proportion" = factor(c("No","Yes")[(deliver_plot_heatmap[,2]>=0.35) + 1], levels=c("Yes","No")),
  "Molecular Class" = factor(c(c("Metabolites","Proteins","Lipids","Oligosaccharides","Peptides","Metabolites"),sapply(rownames(deliver_plot_heatmap)[7:nrow(deliver_plot_heatmap)], molecular_classification)),
                             levels=c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides")),
  col = list(
    "Factor Dominance" = setNames(c(brewer.pal(3,"Set1"),"#CCCCCC"), c("Brand-Region","Region","Time","Residual")),
    "Region Proportion" = setNames(c("yellow","#EEEEEE"), c("Yes","No")),
    "Molecular Class" = setNames(c("#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"), c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides"))
  ),
  annotation_legend_param = list("Factor Dominance" = list(nrow = 2),"Molecular Class" = list(nrow = 2)),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left"
)
rownames(deliver_plot_heatmap)[1:6] = paste0("Total ", rownames(deliver_plot_heatmap)[1:6])
rownames(deliver_plot_heatmap)[7:nrow(deliver_plot_heatmap)] = sapply(rownames(deliver_plot_heatmap)[7:nrow(deliver_plot_heatmap)], map_compound_name)
rownames(deliver_plot_heatmap) = gsub("-[(].*", "", rownames(deliver_plot_heatmap))
rownames(deliver_plot_heatmap) = gsub("Glycosylation-dependent cell adhesion molecule 1", "Protein (UniProt: P80195)", rownames(deliver_plot_heatmap))

svg(paste0(FIGURE_PATH,"\\SuppFig8_RegionDominantRatio_SelectedCompounds.svg"), width=8, height=4)
h = Heatmap(t(deliver_plot_heatmap), name="Proportion of Variance",
            top_annotation = ha,
            column_split = deliver_column_split,
            show_row_names = TRUE, show_column_names = T,
            row_names_gp = gpar(fontsize=10),
            column_names_gp = gpar(fontsize=8),
            column_title_gp = gpar(fontsize=10, fontface="bold"),
            column_names_rot = 75,
            cluster_rows = FALSE, cluster_columns = T,
            cluster_column_slices = F,
            row_names_side = "left",
            show_heatmap_legend = F,
            show_column_dend = F,
            use_raster = T,
            raster_by_magick = TRUE,
            heatmap_legend_param = list(direction = "horizontal"),
            col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")))
draw(h, merge_legend = T, heatmap_legend_side = "bottom",padding = unit(c(3, 5, 3, 3), "mm"))
dev.off()

#All bioactive compounds
bioactivity_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"),header=T,check.names = F)
bioactivity_peptide_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkPeptideBioactivitySummary.csv"),header=T,check.names = F)

all_bioactive_small_molecule = intersect(rownames(tmp_matrix),unique(bioactivity_data$DMD_ID))
all_bioactive_peptide = setdiff(unique(bioactivity_peptide_data$DMD_ID_Peptide),NA)
all_bioactive_protein = unique(bioactivity_peptide_data$DMD_ID_Protein)

deliver_plot_heatmap2 = rbind(
                             tmp_matrix[all_bioactive_small_molecule,],
                             tmp_matrix[all_bioactive_peptide,],
                             tmp_matrix[all_bioactive_protein,])

#Heatmap of the propotion for the whole compounds
ha = HeatmapAnnotation(
  "Factor Dominance" = factor(c("Brand-Region","Region","Time","Residual")[apply(deliver_plot_heatmap2,1,which.max)], levels=c("Brand-Region","Region","Time","Residual")),
  "Region Proportion" = factor(c("No","Yes")[(deliver_plot_heatmap2[,2]>=0.35) + 1], levels=c("Yes","No")),
  "Molecular Class" = factor(sapply(rownames(deliver_plot_heatmap2), molecular_classification),
                             levels=c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides")),
  col = list(
    "Factor Dominance" = setNames(c(brewer.pal(3,"Set1"),"#CCCCCC"), c("Brand-Region","Region","Time","Residual")),
    "Region Proportion" = setNames(c("yellow","#EEEEEE"), c("Yes","No")),
    "Molecular Class" = setNames(c("#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"), c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides"))
  ),
  annotation_legend_param = list("Factor Dominance" = list(nrow = 2),"Molecular Class" = list(nrow = 2)),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left"
)
svg(paste0(FIGURE_PATH,"\\SuppFig24b_HierarchicalVarianceDecomposition_BioactiveCompounds.svg"), width=8, height=2.5)
h = Heatmap(t(deliver_plot_heatmap2), name="Proportion of Variance",
            top_annotation = ha,
            column_split = c("Brand-Region","Region","Time","Residual")[apply(deliver_plot_heatmap2,1,which.max)],
            show_row_names = TRUE, show_column_names = F,
            row_names_gp = gpar(fontsize=10),
            column_names_gp = gpar(fontsize=8),
            column_title_gp = gpar(fontsize=10, fontface="bold"),
            column_names_rot = 75,
            cluster_rows = FALSE, cluster_columns = T,
            cluster_column_slices = F,
            row_names_side = "left",
            show_column_dend = F,
            show_heatmap_legend = F,
            use_raster = T,
            raster_by_magick = TRUE,
            heatmap_legend_param = list(direction = "horizontal"),
            col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")))
draw(h, merge_legend = T, heatmap_legend_side = "bottom",padding = unit(c(3, 5, 3, 3), "mm"))
dev.off()

#How region affect the compound concentrations?
resultB = lapply(colnames(ordered_compound_percentage_nonzero)[union(region_dominate_idx, idx_region_ge_35)], function(x){
  df = data.frame(
    Value = ordered_compound_percentage_nonzero[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand
  )
  
  model = lme4::lmer(Value ~ Region + Time + (1|Brand:Region), data = df)
  emm_region = emmeans(model, "Region")
  # Compare each region to the overall (grand) mean:
  ctr <- contrast(emm_region, method = "eff")  # "effects" contrasts: level - grand mean
  return(ctr)
})

region_diff_emm = t(sapply(resultB, function(x){
  summary(x)$estimate
}))
region_diff_pvalue = t(sapply(resultB, function(x){
  summary(x)$p.value
}))

rownames(region_diff_emm) = colnames(ordered_compound_percentage_nonzero)[union(region_dominate_idx, idx_region_ge_35)]
rownames(region_diff_pvalue) = colnames(ordered_compound_percentage_nonzero)[union(region_dominate_idx, idx_region_ge_35)]
colnames(region_diff_emm) = levels(annotation$Region)
colnames(region_diff_pvalue) = levels(annotation$Region)

region_diff_qvalue = p.adjust(as.vector(region_diff_pvalue), method = "BH")
region_diff_qvalue_matrix = matrix(region_diff_qvalue, nrow=nrow(region_diff_pvalue), ncol=ncol(region_diff_pvalue))
rownames(region_diff_qvalue_matrix) = rownames(region_diff_pvalue)
colnames(region_diff_qvalue_matrix) = colnames(region_diff_pvalue)
                                                                             
idx_row_significant = which(sapply(1:nrow(region_diff_qvalue_matrix), 
                                   function(i){
                                     any(region_diff_qvalue_matrix[i,] < 0.1 & 
                                           abs(region_diff_emm[i,]/mean(ordered_compound_percentage_nonzero[,union(region_dominate_idx, idx_region_ge_35)[i]])) > 0.1)
                                   }
                                   ))

significant_top20_bioactive_compounds = intersect(top20_bioactive_compounds, rownames(region_diff_qvalue_matrix)[idx_row_significant])
significant_bioactive_compounds_only = intersect(bioactive_compounds_only, rownames(region_diff_qvalue_matrix)[idx_row_significant])
significant_top20_compound_only = intersect(top20_compound_only, rownames(region_diff_qvalue_matrix)[idx_row_significant])
significant_other_bioactive_compounds = intersect(setdiff(rownames(deliver_plot_heatmap2), c(significant_top20_bioactive_compounds, significant_bioactive_compounds_only)), rownames(region_diff_qvalue_matrix)[idx_row_significant])


resultB_molecule_class = lapply(c("Oligosaccharides","Peptides"), function(x){
  df = data.frame(
    Value = sum_molecular_class[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand
  )
  
  model = lme4::lmer(Value ~ Region + Time + (1|Brand:Region), data = df)
  emm_region = emmeans(model, "Region")
  # Compare each region to the overall (grand) mean:
  ctr <- contrast(emm_region, method = "eff")  # "effects" contrasts: level - grand mean
  return(ctr)
})

region_diff_emm_molecule_class = t(sapply(resultB_molecule_class, function(x){
  summary(x)$estimate
}))
region_diff_pvalue_molecule_class = t(sapply(resultB_molecule_class, function(x){
  summary(x)$p.value
}))

rownames(region_diff_emm_molecule_class) = c("Oligosaccharides","Peptides")
rownames(region_diff_pvalue_molecule_class) = c("Oligosaccharides","Peptides")
colnames(region_diff_emm_molecule_class) = levels(annotation$Region)
colnames(region_diff_pvalue_molecule_class) = levels(annotation$Region)

region_diff_qvalue_molecule_class = p.adjust(as.vector(region_diff_pvalue_molecule_class), method = "BH")
region_diff_qvalue_matrix_molecule_class = matrix(region_diff_qvalue_molecule_class, nrow=nrow(region_diff_pvalue_molecule_class), ncol=ncol(region_diff_pvalue_molecule_class))
rownames(region_diff_qvalue_matrix_molecule_class) = rownames(region_diff_pvalue_molecule_class)
colnames(region_diff_qvalue_matrix_molecule_class) = colnames(region_diff_pvalue_molecule_class)


#Heatmap of selected compounds
#svg("Figures\\FigS_Heatmap_Selected_Compounds.svg", width = 8, height = 4)
#selected_compounds_top20_bioactive = intersect(union(names(region_dominate_idx), names(idx_region_ge_35)),
#                                               top20_bioactive_compounds)
#selected_compounds_other_bioactive = intersect(union(names(region_dominate_idx), names(idx_region_ge_35)),
#                                               bioactive_compounds_only)
#selected_compounds_other_top20 = intersect(union(names(region_dominate_idx), names(idx_region_ge_35)),
#                                           top20_compound_only)

#matrix_to_plot = t(cbind(sum_molecular_class[,c("Oligosaccharides","Peptides")], 
#                       ordered_compound_percentage_nonzero[,selected_compounds_top20_bioactive],
#                       ordered_compound_percentage_nonzero[,selected_compounds_other_top20],
#                       ordered_compound_percentage_nonzero[,selected_compounds_other_bioactive]
#                       ))
#matrix_to_plot_percentage = t(apply(matrix_to_plot,1,function(x){
#  return((x - mean(x))/mean(x) * 100)
#}))

#rownames(matrix_to_plot_percentage)[1:2] = paste0("Total ", rownames(matrix_to_plot_percentage)[1:2])
#rownames(matrix_to_plot_percentage)[3:nrow(matrix_to_plot_percentage)] = sapply(rownames(matrix_to_plot_percentage)[3:nrow(matrix_to_plot_percentage)], map_compound_name)
#rownames(matrix_to_plot_percentage) = gsub("-[(].*", "", rownames(matrix_to_plot_percentage))
#rownames(matrix_to_plot_percentage) = gsub("Glycosylation-dependent cell adhesion molecule 1", "Protein (UniProt: P80195)", rownames(matrix_to_plot_percentage))

#ha = HeatmapAnnotation(Brand = annotation$Brand,
#                       Time = annotation$Time,
#                       Region = annotation$Region,
#                       
#                       col = list(
#                         Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,desaturate(x,c(0.8)))})),
#                                          levels(annotation$Brand)),
#                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set1"),c("T1","T2","T3")),
#                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
#                       ),
#                       show_legend = F,
#                       annotation_name_side = "left",
#                       annotation_name_rot = c(0,0,0,0),
#                       annotation_name_gp = gpar(fontsize=9),
#                       annotation_legend_param = list(Time = list(nrow = 1),Region = list(nrow = 2)))

#row_split = factor(c(rep("Molecule Class",2), 
#                     rep("Bioactive in top 20", length(selected_compounds_top20_bioactive)), 
#                     rep("Other top 20", length(selected_compounds_other_top20)), 
#                     rep("Other bioactive", length(selected_compounds_other_bioactive))),
#                   levels = c("Molecule Class", "Bioactive in top 20", "Other top 20", "Other bioactive"))

#ra = rowAnnotation(class = row_split,
#                   col = list(class = setNames(c("#1B9E77","#D95F02","#7570B3","#E7298A"), levels(row_split))),
#                   show_legend = F,
#                   show_annotation_name = FALSE)

#h= Heatmap(matrix_to_plot_percentage,
#           row_split = row_split,
#           name = "Difference from Mean (%)",
#           left_annotation = ra,
#           top_annotation = ha,
#           row_names_side="left",
#           show_row_names = T,
#           show_column_names = F,
#           cluster_rows = F,
#           cluster_row_slices = F,
#           cluster_columns = T,
#           show_row_dend = F,
#           show_column_dend = F,
#           row_title = " ",
#           row_gap = unit(1.5, "mm"),
#           row_names_gp = gpar(fontsize=8),
#           row_names_max_width = unit(10, "cm"),
#           #row_title_gp = gpar(fontsize=10, fontface="bold"),
#           #row_title_rot = 0,
#           column_title = "60 Milk Samples",
#           column_title_side = "bottom",
#           column_title_gp = gpar(fontsize=10),
#           col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
#           heatmap_legend_param = list(at=c(-50, -25, 0, 25, 50),labels=c("-50%","-25%","0%","25%","50%"), 
#                                       direction = "horizontal", title_position = "leftcenter")
#)

#draw(h, merge_legend = T, ht_gap=unit(0.5,"cm"), heatmap_legend_side = "bottom")
#dev.off()



svg(paste0(FIGURE_PATH, "\\Fig4a_region_compound_association.svg"), width = 4, height = 5)

matrix_to_plot = t(cbind(sum_molecular_class[,c("Oligosaccharides","Peptides")], 
                         ordered_compound_percentage_nonzero[,selected_compounds_top20_bioactive],
                         ordered_compound_percentage_nonzero[,selected_compounds_other_top20],
                         ordered_compound_percentage_nonzero[,selected_compounds_other_bioactive]
))
matrix_to_plot_percentage = t(apply(matrix_to_plot,1,function(x){
  return((x - mean(x))/mean(x) * 100)
}))
matrix_to_plot_percentage2 = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot_percentage[,col_idx], na.rm=T))
})
significant_ind = matrix(F, nrow=nrow(matrix_to_plot_percentage2), ncol=ncol(matrix_to_plot_percentage2))
rownames(significant_ind) = rownames(matrix_to_plot_percentage2)
colnames(significant_ind) = colnames(matrix_to_plot_percentage2)
significant_ind[intersect(rownames(significant_ind), rownames(region_diff_qvalue_matrix)), 
                colnames(region_diff_qvalue_matrix)] = region_diff_qvalue_matrix[intersect(rownames(significant_ind), rownames(region_diff_qvalue_matrix)), 
                                                                                 colnames(region_diff_qvalue_matrix)] < 0.1
rownames(matrix_to_plot_percentage)[1:2] = paste0("Total ", rownames(matrix_to_plot_percentage)[1:2])
rownames(matrix_to_plot_percentage)[3:nrow(matrix_to_plot_percentage)] = sapply(rownames(matrix_to_plot_percentage)[3:nrow(matrix_to_plot_percentage)], map_compound_name)
rownames(matrix_to_plot_percentage) = gsub("-[(].*", "", rownames(matrix_to_plot_percentage))
rownames(matrix_to_plot_percentage) = gsub("Glycosylation-dependent cell adhesion molecule 1", "Protein (UniProt: P80195)", rownames(matrix_to_plot_percentage))
rownames(matrix_to_plot_percentage2) = rownames(matrix_to_plot_percentage)
rownames(significant_ind) = rownames(matrix_to_plot_percentage)

ha2 = HeatmapAnnotation(Region = factor(regions,levels=regions),
                        col = list(
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F)

ra = rowAnnotation(class = row_split,
                   col = list(class = setNames(c("#1B9E77","#D95F02","#7570B3","#E7298A"), levels(row_split))),
                   show_legend = F,
                   show_annotation_name = FALSE)

h2 = Heatmap(matrix_to_plot_percentage2,
             name = "Difference from Mean (%)",
             row_split = row_split,
             bottom_annotation = ha2,
             left_annotation = ra,
             row_names_side="left",
             show_row_names = T,
             show_column_names = T,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             show_heatmap_legend = F,
             row_title = " ",
             row_names_gp = gpar(fontsize=8),
             row_gap = unit(1.5, "mm"),
             column_title = "",
             column_title_side="bottom",
             column_title_gp = gpar(fontsize=10),
             column_names_gp = gpar(fontsize=9),
             column_names_rot = 75,
             col = colorRamp2(c(-50, -25, 0, 25, 50), c("#0571B0","#92C5DE","#F7F7F7","#F4A582","#CA0020")),
             heatmap_legend_param = list(at=c(-50, -25, 0, 25, 50),labels=c("-50%","-25%","0%","25%","50%")),
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (significant_ind[i, j]) {
                 use_col = "black"
                 if (abs(matrix_to_plot_percentage2[i, j]) > 25) use_col = "white"
                 if (matrix_to_plot_percentage2[i, j] > 0) {
                   grid.text("↑", x, y, 
                             gp = gpar(fontsize = 16, fontface="bold", col=use_col))
                 } else {
                   grid.text("↓", x, y, 
                             gp = gpar(fontsize = 16, fontface="bold", col=use_col))
                 }
               }
             }
)

draw(h2, merge_legend = T, ht_gap=unit(0.5,"cm"), heatmap_legend_side = "bottom")
dev.off()


svg(paste0(FIGURE_PATH,"\\Fig4b_bioactivity_difference_from_mean.svg"), width = 4.5, height = 5)
#Heatmap of bioactivity
plot_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"), header = T, check.names = F)
n_metadata_plot_data = 9
idx_efficacy_plot_data = n_metadata_plot_data + (1:60)
matrix_to_plot = plot_data[,idx_efficacy_plot_data]

selected_bioactivity = c("DPP-IV inhibition - T2D therapy (Nongonierma et al. 2013)",
                         
                         "15-hLO-1 inhibition (AID:887)",
                         
                         "Polymerase Eta inhibition (AID:588591)",
                         "Tdp1 inhibition (AID:485290)",
                         "Polymerase Iota inhibition (AID:720496)",
                         "WRN inhibition (AID:651768)",
                         "15-hLO-2 inhibition (AID:881)",
                         "VDR inhibition - hyperparathyroidism therapy (AID:504847)",
                         "Tau aggregation inhibition - Alzheimer's therapy (AID:1460)",
                         "S. mansoni TGR inhibition (AID:485364)",
                         
                         "RECQ1 inhibition (AID:2549)",
                        
                         "G9a inhibition (AID:504332)",
                         "Malaria parasite viability reduction (AID:504749)"
)
row_split_bioactivity = c("Kappa-casein",
                          "Octadecadienoic acid",
                          rep("Octadecatrienoic acid", 8),
                          "Icosatetraenoic acid",
                          rep("Docosahexaenoic acid", 2))
row_split_bioactivity = factor(row_split_bioactivity, levels = unique(row_split_bioactivity))
rownames(matrix_to_plot) = plot_data$FinalLabel
matrix_to_plot = matrix_to_plot[selected_bioactivity, , drop=F]

bioactivity_class = plot_data$FinalClass[match(rownames(matrix_to_plot), plot_data$FinalLabel)]

ra = rowAnnotation(class = bioactivity_class,
                   col = list(class = setNames(brewer.pal(5,"Set2"), levels(factor(plot_data$FinalClass)))),
                   show_legend = F,
                   show_annotation_name = FALSE)


matrix_to_plot2 = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot[,col_idx], na.rm=T))
})
matrix_to_plot2_bioactivity_percentage = t(apply(matrix_to_plot2,1,function(x){
  return((x - mean(x))/mean(x) * 100)
}))
rownames(matrix_to_plot2_bioactivity_percentage) = gsub("(","\n(", rownames(matrix_to_plot2_bioactivity_percentage), fixed = T)
rownames(matrix_to_plot2_bioactivity_percentage) = gsub("Hypercalcemic hyperparathyroidism \n(AID:504847)", "Hypercalcemic\nhyperparathyroidism\n(AID:504847)", rownames(matrix_to_plot2_bioactivity_percentage), fixed = T)
rownames(matrix_to_plot2_bioactivity_percentage) = gsub("Schistosomiasis, TGR inhibition \n(AID:485364)", "Schistosomiasis,\nTGR inhibition\n(AID:485364)", rownames(matrix_to_plot2_bioactivity_percentage), fixed = T)

h2b = Heatmap(matrix_to_plot2_bioactivity_percentage,
             name = "Difference from Mean (%)",
             row_split = row_split_bioactivity,
             bottom_annotation = ha2,
             right_annotation = ra,
             row_names_side="right",
             show_row_names = T,
             show_column_names = T,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             show_heatmap_legend = T,
             row_title = " ",
             row_names_gp = gpar(fontsize=7.5),
             row_gap = unit(1.5, "mm"),
             column_title = "",
             column_title_side="bottom",
             column_title_gp = gpar(fontsize=10),
             column_names_gp = gpar(fontsize=9),
             column_names_rot = 75,
             col = colorRamp2(c(-20, -10, 0, 10, 20), c("#2B83BA","#ABDDA4","#FFFFBF","#FDAE61","#D7191C")),
             heatmap_legend_param = list(at=c(-20, -10, 0, 10, 20),labels=c("-20%","-10%","0%","10%","20%"))
)
draw(h2b, merge_legend = T, ht_gap=unit(0.5,"cm"), heatmap_legend_side = "bottom")
dev.off()




