#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))

source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

library(ggplot2)
library(RColorBrewer)
library(circlize)

top_10_compound_id = names(colMeans(ordered_compound_percentage)[order(colMeans(ordered_compound_percentage), decreasing = T)][1:10])
df_top_compounds_for_barplot = data.frame(
  percentage = 100*unlist(as.vector(ordered_compound_percentage_wholemilk[,top_10_compound_id])),
  sample_id = c(rep(rownames(ordered_compound_percentage_wholemilk), times = length(top_10_compound_id))),
  Compound = factor(c(sapply(top_10_compound_id, function(x) rep(map_compound_name(x), times = nrow(ordered_compound_percentage_wholemilk)))),
                    levels = c(rev(sapply(top_10_compound_id, function(x) map_compound_name(x))))),
  MolecularClass = factor(c(sapply(top_10_compound_id, function(x) rep(molecular_classification(x), times = nrow(ordered_compound_percentage_wholemilk)))),
                          levels = c("Proteins","Lipids","Metabolites"))
)

df_text = data.frame(
  x = colMeans(ordered_compound_percentage_wholemilk[,top_10_compound_id])*100 + 0.3,
  y = rev(1:10),
  label = sapply(top_10_compound_id, function(x){
    mean_percentage = mean(ordered_compound_percentage_wholemilk[,x])*100
    sd_percentage = sd(ordered_compound_percentage_wholemilk[,x])*100
    return(paste0(map_compound_name(x)," ",round(mean_percentage,2),"% ± ",round(sd_percentage,2),"%"))
  })
)
df_text$x[1] = colMeans(ordered_compound_percentage_wholemilk[,top_10_compound_id[1],drop=F])*30



molecular_class = sapply(colnames(data_pca_tsne), molecular_classification)
sum_top_N_compound_percentage_per_class = lapply(c(10,20,50,100,ncol(ordered_compound_percentage)), function(N){
  top_N_compound_percentage = ordered_compound_percentage[,1:N]
  sum_top_N_compound_percentage_per_class = sapply(unique(molecular_class), function(mol_class){
    mol_ids = colnames(top_N_compound_percentage)[which(molecular_class[colnames(top_N_compound_percentage)] == mol_class)]
    return(rowSums(top_N_compound_percentage[, mol_ids, drop=F]))
  })
  colnames(sum_top_N_compound_percentage_per_class) = unique(molecular_class)
  return(sum_top_N_compound_percentage_per_class)
})

sum_top_N_compound_percentage_per_class_wholemilk = lapply(c(10,20,50,100,ncol(ordered_compound_percentage)), function(N){
  top_N_compound_percentage = ordered_compound_percentage_wholemilk[,1:N]
  sum_top_N_compound_percentage_per_class = sapply(unique(molecular_class), function(mol_class){
    mol_ids = colnames(top_N_compound_percentage)[which(molecular_class[colnames(top_N_compound_percentage)] == mol_class)]
    return(rowSums(top_N_compound_percentage[, mol_ids, drop=F]))
  })
  colnames(sum_top_N_compound_percentage_per_class) = unique(molecular_class)
  return(sum_top_N_compound_percentage_per_class)
})


df_top_compounds_sum_per_class_for_barplot = data.frame(
  percentage = 100*do.call(c,lapply(sum_top_N_compound_percentage_per_class_wholemilk,function(x){as.vector(x)})),
  sample_id = rep(rep(rownames(sum_top_N_compound_percentage_per_class_wholemilk[[1]]), times = ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]])),
                  times=length(sum_top_N_compound_percentage_per_class_wholemilk)),
  Compound = factor(c(rep("Sum of top-10\nCompounds", times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]])*ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]])),
                      rep("Sum of top-20\nCompounds", times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]])*ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]])),
                      rep("Sum of top-50\nCompounds", times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]])*ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]])),
                      rep("Sum of top-100\nCompounds", times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]])*ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]])),
                      rep("Sum of 5,310\nQuantified\nCompounds", times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]])*ncol(sum_top_N_compound_percentage_per_class_wholemilk[[1]]))),
                    levels = c("Sum of 5,310\nQuantified\nCompounds", 
                               "Sum of top-100\nCompounds", 
                               "Sum of top-50\nCompounds",
                               "Sum of top-20\nCompounds",
                               "Sum of top-10\nCompounds")),
  MolecularClass = factor(rep(sapply(colnames(sum_top_N_compound_percentage_per_class_wholemilk[[1]]), function(x) rep(x, times = nrow(sum_top_N_compound_percentage_per_class_wholemilk[[1]]))), 
                              times = length(sum_top_N_compound_percentage_per_class_wholemilk)),
                          levels = c("Proteins","Lipids","Metabolites","Oligosaccharides","Peptides"))
)


df_top_compounds_sum_for_barplot = data.frame(
  percentage = 100*c(as.vector(sum_top_N_compound_percentage_wholemilk)),
  sample_id = c(rep(rownames(sum_top_N_compound_percentage_wholemilk), times = ncol(sum_top_N_compound_percentage_wholemilk))),
  Compound = factor(c(sapply(colnames(sum_top_N_compound_percentage_wholemilk)[1:4], function(x) rep(paste0("Sum of ",gsub("top_","top-",x),"\nCompounds"), times = nrow(sum_top_N_compound_percentage_wholemilk))),
                      rep("Sum of 5,310\nQuantified\nCompounds", times = nrow(sum_top_N_compound_percentage_wholemilk))),
                    levels = c("Sum of 5,310\nQuantified\nCompounds",
                               "Sum of top-100\nCompounds", 
                               "Sum of top-50\nCompounds",
                               "Sum of top-20\nCompounds",
                               "Sum of top-10\nCompounds"))
)


par(mar = c(3, 6, 3, 8), xpd = NA)
svg(paste0(FIGURE_PATH,"\\SuppFig3b.svg"), width = 8, height = 4)
mean_sd = function(x){
  return(c(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x)))
}
label_percentage_of_top_N_compounds = c(sapply(1:5,function(x){
  return(paste0(round(100*(mean(rowSums(sum_top_N_compound_percentage_per_class_wholemilk[[x]]))),1),"%"))
}))
label_percentage_of_top_N_compounds = paste0(label_percentage_of_top_N_compounds, " (",
                                             c(sapply(1:5,function(x){
                                               return(paste0(round(100*(mean(rowSums(sum_top_N_compound_percentage_per_class_wholemilk[[x]]))/mean(sum_top_N_compound_percentage_wholemilk[,5])),1),"%"))
                                             })),")")



saccharide_label = paste0("",round(sapply(1:5, function(i){colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Metabolites"]*100}),1),"%")
saccharide_pos = sapply(1:5, function(i){
  sum(colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])[c("Peptides","Oligosaccharides")])*100 +
    colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Metabolites"]*100/2
})

fat_label = paste0("",round(sapply(1:5, function(i){colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Lipids"]*100}),1),"%")
fat_pos = sapply(1:5, function(i){
  sum(colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])[c("Peptides","Oligosaccharides","Metabolites")])*100 +
    colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Lipids"]*100/2
})

protein_label = paste0("",round(sapply(1:5, function(i){colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Proteins"]*100}),1),"%")
protein_pos = sapply(1:5, function(i){
  sum(colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])[c("Peptides","Oligosaccharides","Metabolites","Lipids")])*100 +
    colMeans(sum_top_N_compound_percentage_per_class_wholemilk[[i]])["Proteins"]*100/2
})

df_text = data.frame(
  x = c(rep(13.55,5), saccharide_pos, fat_pos, protein_pos),
  y = rep(c(5,4,3,2,1),4),
  label = c(label_percentage_of_top_N_compounds,
            saccharide_label,
            fat_label,
            protein_label)
)
p = ggplot(df_top_compounds_sum_per_class_for_barplot, aes(y=Compound, x=percentage)) +
  geom_bar(stat="summary", fun = "mean", position = "stack", aes(fill=MolecularClass), width=0.5) +
  #stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.3, position = position_stack(vjust = 0.5)) + # For error bars (mean +/- standard error)
  #geom_jitter(width = 0, height = 0.4, alpha = 0.1, size=0.2, show.legend = F) +
  theme_classic() +
  labs(title="", y="", x="% Content (w/w)", fill="") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right", 
        legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=8))+
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"),guide = guide_legend(ncol = 1)) +
  scale_x_continuous(limits = c(0,17), breaks=c(0,5,10,15), expand = c(0,0), 
                     name = "% Content (w/w)",
                     sec.axis = sec_axis(
                       # Inverse transformation for the secondary axis (right)
                       ~. * 100 / 11.5, 
                       name = "% Content of Quantified Compounds (w/w)", breaks=c(0,50,100)
                     )) + 
  geom_jitter(data=df_top_compounds_sum_for_barplot, aes(y=Compound, x=percentage), fill="black", width = 0, height = 0.3, alpha = 0.1, size=0.5) + 
  #geom_bar(data=df_top_compounds_sum_for_barplot, aes(y=Compound, x=percentage), fill=NA, color="black",
  #         linetype="dashed", linewidth=0.5, stat="summary", fun = "mean") +
  geom_errorbar(data=df_top_compounds_sum_for_barplot, aes(y=Compound, x=percentage),
                stat="summary", fun.data = mean_sd, width = 0.3, color="black",
                linewidth=0.5) + 
  geom_text(data=df_text, aes(x=x, y=y, label=label), size=2.8) 
print(p)
dev.off()




#Prepare sankey input
#(Current class, Current Node, Next Class, Next Node, Value)
percentage_class_all = sapply(levels(df_top_compounds_sum_per_class_for_barplot$MolecularClass), function(x){
  df_sub = df_top_compounds_sum_per_class_for_barplot[which(df_top_compounds_sum_per_class_for_barplot$MolecularClass == x), ]
  df_sub = df_sub$percentage[which(df_sub$Compound == "Sum of 5,310\nQuantified\nCompounds")]
  return(mean(df_sub))
})
percentage_water = 100 - sum(percentage_class_all)

pencentage_top_10_compound = sapply(levels(df_top_compounds_for_barplot$Compound), function(x){
  df_sub = df_top_compounds_for_barplot[which(df_top_compounds_for_barplot$Compound == x), ]
  return(mean(df_sub$percentage))
})
pencentage_top_10_compound = pencentage_top_10_compound[order(pencentage_top_10_compound, decreasing = T)]

pencentage_other_protein = percentage_class_all["Proteins"] - sum(pencentage_top_10_compound[c(2,3,4,6,8,10)])
pencentage_other_lipid = percentage_class_all["Lipids"] - sum(pencentage_top_10_compound[c(5,7,9)])
pencentage_other_metabolite = percentage_class_all["Metabolites"] - sum(pencentage_top_10_compound[c(1)])


library('plotrix')
#Need to plot from outside to inside based on this mechanism
#0 circle: blank
iniR = 0.2

pencentage_other_protein = percentage_class_all["Proteins"] - sum(pencentage_top_10_compound[c(2,3,4,6,8,10)])
pencentage_other_lipid = percentage_class_all["Lipids"] - sum(pencentage_top_10_compound[c(5,7,9)])
pencentage_other_metabolite = percentage_class_all["Metabolites"] - sum(pencentage_top_10_compound[c(1)])

#Base
svg(paste0(FIGURE_PATH,"\\Fig2a_piechart_base.svg"), width = 6, height = 6)
pie(1, radius=iniR, col=c('white'), border = NA, labels='')
start_pt = 275/360*2*pi
#Second Layer
floating.pie(0,0,c(c(pencentage_top_10_compound[c(2,3,4,6,8,10)],
                     pencentage_other_protein),
                   c(pencentage_top_10_compound[c(5,7,9)],
                     pencentage_other_lipid),
                   c(pencentage_top_10_compound[c(1)],
                     pencentage_other_metabolite),
                   percentage_class_all[c(4,5)]
),radius=3.5*iniR, col=c( c("#87BAB1","#63A79C","#3D9487","#277F73","#086A5F","#06554C","#AAEFE3"),
                        c("#E2E2A4","#C8C866","#AEAE00","#FFFFD8"),
                        c("#A9A6BE","#DEDAFB"),
                        rep("white",2)),border=NA,start=start_pt)
#First Layer
floating.pie(0,0,c(percentage_class_all[c(1,2,3,4,5)]),radius=3*iniR, col=c("#8DD3C7","#FFFFB3","#BEBADA","#DE9AEA","#FDB462"),border=NA,start=start_pt)
#Inner Circle
#floating.pie(0,0, 1, startpos=0, radius=2*iniR, col="white", border = NA)
dev.off()




library(ComplexHeatmap)
library(colorspace)


#Correlation of samples
correlation_matrix_top10 = cor(t(ordered_compound_percentage[,1:10]), method = "pearson")
order_cor = hclust(as.dist(1-correlation_matrix_top10))$order
correlation_matrix_top10 = correlation_matrix_top10[order_cor, order_cor]

library(ComplexHeatmap)
annotation = data.frame(
  Region = product_data$`Region Code`[match(rownames(correlation_matrix_top10), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(rownames(correlation_matrix_top10), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(rownames(correlation_matrix_top10), product_data$`DMD ID`)]
)
rownames(annotation) = rownames(annotation)
annotation$Region = factor(regions[annotation$Region],levels=regions)
annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))
annotation$Brand = factor(c("Brand A","Brand B")[annotation$Brand], levels=c("Brand A","Brand B"))

svg(paste0(FIGURE_PATH,"\\SuppFig16a_correlation_top10.svg"), width = 3.5, height = 3.5)

ha = HeatmapAnnotation(Brand = annotation$Brand,
                       Time = annotation$Time,
                       Region = annotation$Region,
                       col = list(
                         Brand = setNames(c("#E41A1C","#377EB8"),c("Brand A","Brand B")),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set2"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize=9))

h = Heatmap(correlation_matrix_top10,
            name = "Correlation",
            bottom_annotation = ha,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            cluster_columns = F,
            show_row_dend = F,
            show_column_dend = F,
            row_title = "",
            column_title = "",
            column_title_gp = gpar(fontsize=8),
            col = colorRamp2(c(0.8, 0.85, 0.90, 0.95, 1), c("#253494","#2C7FB8","#41B6C4","#A1DAB4","#FFFFCC")),
            heatmap_legend_param = list(at=c(0.8, 0.85, 0.90, 0.95, 1),labels=c("0.8","0.85","0.9","0.95","1")),
            rect_gp = gpar(type="none"),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (i >= j){
                grid.rect(x, y, width, height, gp = gpar(fill = fill, col = NA))
              }
            }
)
draw(h, heatmap_legend_side = "right")
dev.off()



correlation_matrix_all = cor(t(ordered_compound_percentage), method = "pearson")
order_cor = hclust(as.dist(1-correlation_matrix_all))$order
correlation_matrix_all = correlation_matrix_all[order_cor, order_cor]

library(ComplexHeatmap)
annotation = data.frame(
  Region = product_data$`Region Code`[match(rownames(correlation_matrix_all), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(rownames(correlation_matrix_all), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(rownames(correlation_matrix_all), product_data$`DMD ID`)]
)
rownames(annotation) = rownames(annotation)
annotation$Region = factor(regions[annotation$Region],levels=regions)
annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))
annotation$Brand = factor(c("Brand A","Brand B")[annotation$Brand], levels=c("Brand A","Brand B"))


svg(paste0(FIGURE_PATH,"\\SuppFig16b_Correlation_AllCompounds.svg"), width = 3.5, height = 3.5)
ha = HeatmapAnnotation(Brand = annotation$Brand,
                       Time = annotation$Time,
                       Region = annotation$Region,
                       col = list(
                         Brand = setNames(c("#E41A1C","#377EB8"),c("Brand A","Brand B")),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set2"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize=3))

h = Heatmap(correlation_matrix_all,
            name = "Correlation",
            bottom_annotation = ha,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            cluster_columns = F,
            show_row_dend = F,
            show_column_dend = F,
            row_title = "",
            column_title = "",
            column_title_gp = gpar(fontsize=8),
            col = colorRamp2(c(0.8, 0.85, 0.90, 0.95, 1), c("#253494","#2C7FB8","#41B6C4","#A1DAB4","#FFFFCC")),
            heatmap_legend_param = list(at=c(0.8, 0.85, 0.90, 0.95, 1),labels=c("0.8","0.85","0.9","0.95","1")),
            rect_gp = gpar(type="none"),
            show_heatmap_legend = F,
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (i >= j){
                grid.rect(x, y, width, height, gp = gpar(fill = fill, col = NA))
              }
            }
)
draw(h, heatmap_legend_side = "right")
dev.off()





#Stacked bar for the top 10 compounds
df_stacked_bar_top10 = data.frame(
  percentage = c(100*as.vector(as.matrix(ordered_compound_percentage_wholemilk[,top_10_compound_id])),
                 100*rowSums(ordered_compound_percentage_wholemilk[,setdiff(colnames(ordered_compound_percentage_wholemilk), top_10_compound_id)])),
  sample_id = c(rep(rownames(ordered_compound_percentage_wholemilk), times = length(top_10_compound_id) + 1)),
  Compound = factor(c(sapply(top_10_compound_id, function(x) rep(map_compound_name(x), times = nrow(ordered_compound_percentage_wholemilk))),
                      rep("Other Compounds", times = nrow(ordered_compound_percentage_wholemilk))),
                    levels = c("Other Compounds",rev(sapply(top_10_compound_id, function(x) map_compound_name(x)))[c(1,3,6:9,2,4,5,10)]))
)

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


matrix_stacked_bar_annotation = cbind(ordered_compound_percentage[,top_10_compound_id], rowSums(ordered_compound_percentage[,setdiff(colnames(ordered_compound_percentage), top_10_compound_id)]))
matrix_stacked_bar_annotation = matrix_stacked_bar_annotation[,(c(1,6,7,9,2:5,8,10,11))]
#ComplexHeatmap without heatmap body -- just add the annotation with stacked bar plot here.
ha = HeatmapAnnotation(
                       "Content (% w/w)" = anno_barplot(100*matrix_stacked_bar_annotation, height = unit(6, "cm"), border=F, 
                                        gp = gpar(fill= c(rev(c("#87BAB1","#63A79C","#3D9487","#277F73","#086A5F","#06554C","#E2E2A4","#C8C866","#AEAE00","#A9A6BE")),"grey"), col=NA),
                                        axis_param = list(at = seq(0,14,2)),ylim=c(0,14)),
                       Brand = annotation$Brand,
                       Time = annotation$Time,
                       Region = annotation$Region,
                       col = list(
                         Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,darken(x,c(0.1)))})),
                                          levels(annotation$Brand)),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set1"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_rot = c(90,0,0,0),
                       annotation_name_gp = gpar(fontsize=9))

h = Heatmap(matrix(NA, nrow=0, ncol=nrow(ordered_compound_percentage)), # Empty matrix for the heatmap body
            name = "Test",
            top_annotation = ha,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            cluster_columns = F,
            show_row_dend = F,
            show_column_dend = F,
            row_title = "",
            column_title = "",
            column_title_gp = gpar(fontsize=8),
            rect_gp = gpar(type="none"),
            show_heatmap_legend = F
)

svg(paste0(FIGURE_PATH,"\\SuppFig3a_stackedbar_top10compounds.svg"), width = 6, height = 4)
draw(h, heatmap_legend_side = "right")
dev.off()


#Violin map for the top 10 compounds and others (Normalized by their means)
matrix_normalized = apply(matrix_stacked_bar_annotation, 2, function(x){
  return((x-mean(x))/mean(x))
})


library(patchwork)
df_boxplot_top10 = data.frame(
  content = 100*as.vector(as.matrix(matrix_stacked_bar_annotation)),
  diff_to_mean = 100*as.vector(as.matrix(matrix_normalized)),
  sample_id = c(rep(rownames(matrix_normalized), times = ncol(matrix_normalized))),
  Compound = factor(c(sapply(colnames(matrix_normalized)[1:10], function(x) rep(map_compound_name(x), times = nrow(matrix_normalized))),
                      rep("Other Compounds", times = nrow(matrix_normalized))),
                    levels = c("Other Compounds",rev(sapply(colnames(matrix_normalized)[1:10], function(x) map_compound_name(x)))))
)
#levels(df_boxplot_top10$Compound) = (levels(df_boxplot_top10$Compound) )
plot_color = (c("grey","#87BAB1","#63A79C","#3D9487","#277F73","#086A5F","#06554C","#E2E2A4","#C8C866","#AEAE00","#A9A6BE"))

mean_sd = function(x){
  return(c(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x)))
}

h_bar = ggplot(df_boxplot_top10, aes(y=Compound, x=content, fill=Compound)) +
  geom_bar(stat="summary", fun = "mean", width=0.5) +
  theme_classic() +
  labs(title="", x="% Content (w/w)", y="", fill="Compound") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text.x = element_text(hjust = 1),
        axis.title.x = element_text(size=10), axis.text = element_text(size=10)) +
  scale_fill_manual(values = plot_color) + 
  scale_x_continuous(limits = c(0, 6), breaks=seq(0,6,1), expand = c(0,0)) + 
  geom_jitter(width = 0, height = 0.3, alpha = 0.1, size=0.5, show.legend = F) + 
  geom_errorbar(stat="summary", fun.data = mean_sd, width = 0.3, color="black", linewidth=0.5) + 
  geom_text(data=df_boxplot_top10 %>% group_by(Compound) %>% summarise(mean_content = mean(content)), 
            aes(x=mean_content + 0.8, y=Compound, label=paste0(round(mean_content,2),"%")), size=2.8)
  

h_box = ggplot(df_boxplot_top10, aes(y=Compound, x=diff_to_mean, fill=Compound)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic() +
  labs(title="", x="% Relative Difference\nto Mean", y="", fill="Compound") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text.x = element_text(hjust = 1)) +
  scale_fill_manual(values = plot_color) + 
  scale_x_continuous(limits = c(-60, 60)) + 
  geom_jitter(width = 0, height = 0.3, alpha = 0.2, size=0.8, show.legend = F) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=10))

#Variance bar
rsd_func = function(x){
  return(sd(x)/mean(x))
}
h_bar_rsd = ggplot(df_boxplot_top10, aes(y=Compound, x=content, fill=Compound)) +
  geom_bar(stat="summary", fun = rsd_func, width=0.5) +
  theme_classic() +
  labs(title="", x="Relative\nStandard Deviation", y="", fill="Compound") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text.x = element_text(hjust = 1)) +
  scale_fill_manual(values = plot_color) + 
  scale_x_continuous(limits = c(0, 0.5), breaks=seq(0,0.5,0.25), oob = scales::oob_keep, expand = c(0,0)) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=10)) + 
  geom_text(data=df_boxplot_top10 %>% group_by(Compound) %>% summarise(rsd = rsd_func(content)), 
            aes(x=rsd + 0.05, y=Compound, label=paste0(round(rsd,2))), size=2.8)

h = h_bar + h_box + h_bar_rsd + plot_layout(ncol = 3, widths = c(1,1,0.5)) 

svg(paste0(FIGURE_PATH,"\\SuppFig4.svg"), width = 8, height = 4)
print(h)
dev.off()




library(dplyr)
library(car)

#The significant test results for determining significantly higher/lower variance of each compound compared to others is added to the Supp. Fig 4c manually.
leveneTest(content ~ Compound, data = df_boxplot_top10)

significant_results = lapply(levels(df_boxplot_top10$Compound), function(target){
  other_compounds <- setdiff(unique(df_boxplot_top10$Compound), target)
  
  pairwise_results <- lapply(other_compounds, function(comp) {
    
    tmp <- df_boxplot_top10 %>%
      filter(Compound %in% c(target, comp))
    
    test <- leveneTest(content ~ Compound, data = tmp)
    
    data.frame(
      target = target,
      comparison = comp,
      F_value = test$`F value`[1],
      p_value = test$`Pr(>F)`[1]
    )
  })
  
  pairwise_results <- bind_rows(pairwise_results) %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      significant = p_adjusted < 0.05
    )
  
  return(pairwise_results)
})
