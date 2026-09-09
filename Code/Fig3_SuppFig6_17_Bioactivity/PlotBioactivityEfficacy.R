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
library(ComplexHeatmap)
library(scales)
library(colorspace)

if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"))){
  stop("Please run the script 'MergeSmallMoleculeAndPeptideBioactivity.r' with necessary manual labeling first to prepare the bioactivity efficacy data for plotting.")
}

#Heatmap
plot_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary_ManualLabel.csv"), header = T, check.names = F)
n_metadata_plot_data = 9
idx_efficacy_plot_data = n_metadata_plot_data + (1:60)
idx_invivo_efficacy_plot_data = n_metadata_plot_data + (61:120)
#a. In vitro

annotation = data.frame(
  Region = product_data$`Region Code`[match(colnames(plot_data[,idx_efficacy_plot_data]), product_data$`DMD ID`)],
  Time = product_data$`Time Point Code`[match(colnames(plot_data[,idx_efficacy_plot_data]), product_data$`DMD ID`)],
  Brand = product_data$`Brand Code`[match(colnames(plot_data[,idx_efficacy_plot_data]), product_data$`DMD ID`)]
)
rownames(annotation) = rownames(annotation)
annotation$Brand = factor(paste0("Brand ",annotation$Region, c("A","B")[annotation$Brand]), 
                          levels=paste0("Brand ",sapply(1:10,function(x){rep(x,2)}), rep(c("A","B"),10)))
annotation$Region = factor(regions[annotation$Region],levels=regions)
annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))


row_split = factor(plot_data$FinalClass, levels = unique(plot_data$FinalClass))
levels(row_split) = c("Anticancer","Disease\nTherapy","Anti-\noxidation","Anti-\nmicrobial","Antivirus")

svg(paste0(FIGURE_PATH,"\\SuppFig6a_InvitroBioactivity.svg"), width = 8, height = 5)
#ramp <- colour_ramp(brewer.pal(9, "Oranges"))

matrix_to_plot = plot_data[,idx_efficacy_plot_data]
rownames(matrix_to_plot) = plot_data$FinalLabel
ha = HeatmapAnnotation(Brand = annotation$Brand,
                       Time = annotation$Time,
                       Region = annotation$Region,
                       "Average Efficacy (%)" = anno_barplot(colMeans(matrix_to_plot), height = unit(1, "cm"), border=F, 
                                                         gp = gpar(fill="grey", col=NA),
                                                         axis_param = list(at = c(0, 50, 80)),ylim=c(0,85)),
                       col = list(
                         Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,darken(x,c(0.5)))})),
                                                  levels(annotation$Brand)),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set1"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_rot = c(0,0,0,0),
                       annotation_name_gp = gpar(fontsize=9))

h= Heatmap(matrix_to_plot,
           row_split = row_split,
           name = "In vitro Efficacy (%)",
           top_annotation = ha,
           row_names_side="left",
           show_row_names = T,
           show_column_names = F,
           cluster_rows = F,
           cluster_row_slices = F,
           cluster_columns = T,
           show_row_dend = F,
           show_column_dend = F,
           row_title = " ",
           row_gap = unit(1.5, "mm"),
           row_names_gp = gpar(fontsize=7),
           row_names_max_width = unit(10, "cm"),
           #row_title_gp = gpar(fontsize=10, fontface="bold"),
           #row_title_rot = 0,
           column_title = "60 Milk Samples",
           column_title_side = "bottom",
           column_title_gp = gpar(fontsize=10),
           col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
           heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%"), 
                                       direction = "horizontal", title_position = "leftcenter")
)

matrix_to_plot2 = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot[,col_idx], na.rm=T))
})

ha2 = HeatmapAnnotation(Dummy1 = factor(regions,levels=regions),
                        Dummy2 = factor(regions,levels=regions),
                        Region = factor(regions,levels=regions),
                        "Average Efficacy" = anno_barplot(colMeans(matrix_to_plot2), height = unit(1, "cm"), border=F, axis=F,
                                                          gp = gpar(fill="grey", col=NA)),
                        col = list(
                          Dummy1 = setNames(rep("white",length(regions)),regions),
                          Dummy2 =  setNames(rep("white",length(regions)),regions),
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F)

h2 = Heatmap(matrix_to_plot2,
             name = "In vitro Efficacy (%)",
             row_split = row_split,
             top_annotation = ha2,
             row_names_side="left",
             show_row_names = F,
             show_column_names = F,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             row_title = " ",
             row_gap = unit(1.5, "mm"),
             column_title = "Region Avg.",
             column_title_side="bottom",
             column_title_gp = gpar(fontsize=10),
             col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%"))
)

draw(h+h2, merge_legend = T, ht_gap=unit(0.2,"cm"), heatmap_legend_side = "bottom")
dev.off()



svg(paste0(FIGURE_PATH,"\\SuppFig6b_InvivoBioactivity.svg"), width = 8, height = 5)

matrix_to_plot = plot_data[,idx_invivo_efficacy_plot_data]
rownames(matrix_to_plot) = plot_data$FinalLabel
ha = HeatmapAnnotation(Brand = annotation$Brand,
                       Time = annotation$Time,
                       Region = annotation$Region,
                       "Average Efficacy (%)" = anno_barplot(round(colMeans(matrix_to_plot),1), height = unit(1, "cm"), border=F, 
                                                         gp = gpar(fill="grey", col=NA),
                                                         axis_param = list(at = c(0, 50, 80)),
                                                         add_numbers = T,
                                                         numbers_gp = gpar(fontsize = 6.5, col = "black"),
                                                         numbers_rot=90,
                                                         numbers_offset = unit(0.5,"mm"),ylim=c(0,50)),
                       col = list(
                         Brand = setNames(as.vector(sapply(brewer.pal(length(unique(annotation$Region)),"Set3"), function(x){c(x,darken(x,c(0.5)))})),
                                           levels(annotation$Brand)),
                         Time = setNames(brewer.pal(length(unique(annotation$Time)),"Set1"),c("T1","T2","T3")),
                         Region = setNames(brewer.pal(length(unique(annotation$Region)),"Set3"),regions)
                       ),
                       show_legend = F,
                       annotation_name_side = "left",
                       annotation_name_rot = c(0,0,0,0),
                       annotation_name_gp = gpar(fontsize=9))


h= Heatmap(matrix_to_plot,
           row_split = row_split,
           name = "In vivo Efficacy (%)",
           top_annotation = ha,
           row_names_side="left",
           show_row_names = T,
           show_column_names = F,
           cluster_rows = F,
           cluster_row_slices = F,
           cluster_columns = T,
           show_row_dend = F,
           show_column_dend = F,
           row_title = " ",
           row_gap = unit(1.5, "mm"),
           row_names_gp = gpar(fontsize=7),
           row_names_max_width = unit(10, "cm"),
           #row_title_gp = gpar(fontsize=10, fontface="bold"),
           #row_title_rot = 0,
           column_title = "60 Milk Samples",
           column_title_side = "bottom",
           column_title_gp = gpar(fontsize=10),
           col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
           heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%"), 
                                       direction = "horizontal", title_position = "leftcenter")
)

matrix_to_plot2 = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot[,col_idx], na.rm=T))
})
ha2 = HeatmapAnnotation(Dummy1 = factor(regions,levels=regions),
                        Dummy2 = factor(regions,levels=regions),
                        Region = factor(regions,levels=regions),
                        "Average Efficacy" = anno_barplot(round(colMeans(matrix_to_plot2),1), height = unit(1, "cm"), border=F, axis=F,
                                                          gp = gpar(fill="grey", col=NA),
                                                          axis_param = list(at = c(0, 50, 80)),ylim=c(0,50),
                                                          add_numbers = T,
                                                          numbers_gp = gpar(fontsize = 6.5, col = "black"),
                                                          numbers_rot=90,
                                                          numbers_offset = unit(0.5,"mm")),
                        col = list(
                          Dummy1 = setNames(rep("white",length(regions)),regions),
                          Dummy2 =  setNames(rep("white",length(regions)),regions),
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F)

h2 = Heatmap(matrix_to_plot2,
             name = "In vivo Efficacy (%)",
             row_split = row_split,
             top_annotation = ha2,
             row_names_side="left",
             show_row_names = F,
             show_column_names = F,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             row_title = " ",
             row_gap = unit(1.5, "mm"),
             column_title = "Region Avg.",
             column_title_side="bottom",
             column_title_gp = gpar(fontsize=10),
             col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%"))
)

draw(h+h2, merge_legend = T, ht_gap=unit(0.2,"cm"), heatmap_legend_side = "bottom")
dev.off()



svg(paste0(FIGURE_PATH,"\\Fig3a_RegionBioactivity.svg"), width = 7, height = 5.2)

matrix_to_plot = plot_data[,idx_efficacy_plot_data]
rownames(matrix_to_plot) = plot_data$FinalLabel
matrix_to_plot2a = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot[,col_idx], na.rm=T))
})
ha2 = HeatmapAnnotation(Region = factor(regions,levels=regions),
                        "Average Efficacy (%)" = anno_barplot(colMeans(matrix_to_plot2a), height = unit(1, "cm"), border=F, axis=T,
                                                          gp = gpar(fill="grey", col=NA),
                                                          axis_param = list(at = c(0, 50, 80)), ylim=c(0,85)),
                        col = list(
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "left",
                        annotation_name_gp = gpar(fontsize=9),
                        annotation_name_rot=0,
                        show_annotation_name = T,
                        show_legend = F)


in_vitro_summary_annotation =  rowAnnotation(in_vitro_avg = anno_text(sprintf("%.2f",rowMeans(matrix_to_plot2a)), 
                                                                      location = 0.5, 
                                                                      just = "center",
                                                                      gp = gpar(fill = rep(2:4, each = 4), col = "white", border = "black"),
                                                                      width = max_text_width(month.name)*1.2))


h2a = Heatmap(matrix_to_plot2a,
             name = "In vitro Efficacy (%)",
             row_split = row_split,
             top_annotation = ha2,
             row_names_side="left",
             show_row_names = T,
             show_column_names = F,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             row_names_gp = gpar(fontsize=7),
             row_title = " ",
             row_gap = unit(1.5, "mm"),
             column_title = "Estimated\nin vitro efficacy (%)",
             column_title_side="top",
             column_title_gp = gpar(fontsize=10, fontface = "bold"),
             col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%")),
             show_heatmap_legend = F
)
avg_invitro = matrix(rowMeans(matrix_to_plot2a), ncol=1)
h2a_addl = Heatmap(avg_invitro,
                   name = "Avg. In vitro Efficacy (%)",
                   row_split = row_split,
                   show_row_names = F,
                   show_column_names = F,
                   cluster_row_slices = F,
                   cluster_rows = F,
                   cluster_columns = T,
                   show_row_dend = F,
                   show_column_dend = F,
                   row_title = " ",
                   row_gap = unit(1.5, "mm"),
                   column_title = "",
                   column_title_gp = gpar(fontsize=8),
                   col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
                   heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%")),
                   show_heatmap_legend = F,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (avg_invitro[i, j] > 50){
                       grid.text(sprintf("%.1f", avg_invitro[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#DDDDDD"))
                     }else{
                       grid.text(sprintf("%.1f", avg_invitro[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#333333"))
                     }
                   },
                   width=unit(0.8,"cm")
                   )

matrix_to_plot = plot_data[,idx_invivo_efficacy_plot_data]
rownames(matrix_to_plot) = plot_data$FinalLabel
matrix_to_plot2b = sapply(levels(annotation$Region), function(region){
  col_idx = which(annotation$Region == region)
  return(rowMeans(matrix_to_plot[,col_idx], na.rm=T))
})
ha2b = HeatmapAnnotation(Region = factor(regions,levels=regions),
                        "Average Efficacy" = anno_barplot(round(colMeans(matrix_to_plot2b),1), height = unit(1, "cm"), border=F, axis=F,
                                                          gp = gpar(fill="grey", col=NA),
                                                          axis_param = list(at = c(0, 50, 80)), ylim=c(0,85),
                                                          add_numbers = T,
                                                          numbers_gp = gpar(fontsize = 7, col = "black"),
                                                          numbers_rot=60),
                        col = list(
                          Region = setNames(brewer.pal(length(regions),"Set3"),regions)
                        ),
                        annotation_name_side = "right",
                        show_annotation_name = F,
                        show_legend = F)

h2b = Heatmap(matrix_to_plot2b,
             name = "In vivo Efficacy (%)",
             row_split = row_split,
             top_annotation = ha2b,
             row_names_side="left",
             show_row_names = F,
             show_column_names = F,
             cluster_row_slices = F,
             cluster_rows = F,
             cluster_columns = T,
             show_row_dend = F,
             show_column_dend = F,
             row_title = " ",
             row_gap = unit(1.5, "mm"),
             column_title = "Estimated\nin vivo efficacy (%)",
             column_title_side="top",
             column_title_gp = gpar(fontsize=10, fontface = "bold"),
             col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
             heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%")),
             show_heatmap_legend = F
)

avg_invivo = matrix(rowMeans(matrix_to_plot2b), ncol=1)
rownames(avg_invivo) = rownames(matrix_to_plot2b)
h2b_addl = Heatmap(avg_invivo,
                   name = "Avg. In vivo Efficacy (%)",
                   row_split = row_split,
                   show_row_names = F,
                   show_column_names = F,
                   cluster_row_slices = F,
                   cluster_rows = F,
                   cluster_columns = T,
                   show_row_dend = F,
                   show_column_dend = F,
                   row_title = " ",
                   row_gap = unit(1.5, "mm"),
                   column_title = "",
                   column_title_gp = gpar(fontsize=8),
                   col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
                   heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%")),
                   show_heatmap_legend = F,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (avg_invivo[i, j] > 50){
                       grid.text(sprintf("%.1f", avg_invivo[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#DDDDDD"))
                     }else{
                       grid.text(sprintf("%.1f", avg_invivo[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#333333"))
                     }
                   },
                   width=unit(0.8,"cm")
)

avg_invivo_extend = plot_data[,paste0("Avg_Invivo_",c(0.01,0.1,10,100),"x")]
colnames(avg_invivo_extend) = paste0(c(0.01,0.1,10,100),"X")
rownames(avg_invivo_extend) = rownames(matrix_to_plot2b)
h2b_addl_extend = Heatmap(avg_invivo_extend,
                          name = "Avg. In vivo Efficacy (%)",
                          row_split = row_split,
                          show_row_names = F,
                          show_column_names = F,
                          cluster_row_slices = F,
                          cluster_rows = F,
                          cluster_columns = F,
                          show_row_dend = F,
                          show_column_dend = F,
                          row_title = " ",
                          row_gap = unit(1.5, "mm"),
                          top_annotation = HeatmapAnnotation(
                            foo=anno_text(colnames(avg_invivo_extend),gp=gpar(fontsize=7),rot=0,location=0.5,just="center")
                            ),
                          column_title = "Estimated in vivo efficacy (%)\n(In-serum concentration variation)",
                          column_title_gp = gpar(fontsize=6.5, fontface = "bold"),
                          col = colorRamp2(c(0, 25, 50, 75, 100), c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026")),
                          heatmap_legend_param = list(at=c(0, 25, 50, 75, 100),labels=c("0%","25%","50%","75%","100%")),
                          show_heatmap_legend = F,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if (avg_invivo_extend[i, j] > 50){
                              grid.text(sprintf("%.1f", avg_invivo_extend[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#DDDDDD"))
                            }else{
                              grid.text(sprintf("%.1f", avg_invivo_extend[i, j]), x, y, gp = gpar(fontsize = 8, fontface="bold", col="#333333"))
                            }
                          },
                          width=unit(3.6,"cm")
)

draw(h2a+h2a_addl+h2b+h2b_addl+h2b_addl_extend, merge_legend = T, ht_gap=unit(c(0.05,0.1,0.05,0.15),"cm"), heatmap_legend_side = "bottom")
dev.off()


#BipartiteGraph (Note: Only plotting the graph itself. The text labels are added manually)
final_label_with_newline = plot_data$FinalLabel
to_label_modified = sapply(1:nrow(plot_data),function(i){
  if (plot_data$CompoundClass[i] == "Peptides"){
    return(paste0("Peptide: ", plot_data$Compound[i]))
  }
  return(gsub("-[(].*","", plot_data$Compound[i]))
})

df_chord_invitro = data.frame(from=to_label_modified, to=final_label_with_newline, value=rowMeans(plot_data[,idx_efficacy_plot_data]))

#grid.col <- setNames(rainbow(length(c(df_chord_invitro$from, df_chord_invitro$to))), unique(c(df_chord_invitro$from, df_chord_invitro$to)))
compound_class_color = c("#8DD3C7","#FFFFB3","#BEBADA","#FDB462")
names(compound_class_color) = c("Proteins","Lipids","Metabolites","Peptides")
grid.colA = brewer.pal(6,"Set2")[as.numeric(factor(plot_data$FinalClass))]
names(grid.colA) = final_label_with_newline
grid.colB = compound_class_color[as.character(plot_data$CompoundClass)]
names(grid.colB) = to_label_modified
grid.col = c(grid.colA, grid.colB)

svg(paste0(FIGURE_PATH,"\\Fig3b_BipartiteGraph.svg"), width = 10, height = 20)
#par(mar = c(3, 3, 3, 3), mfrow = c(1, 1), xpd = NA, pty="m")
circos.par(gap.after = c(rep(0.05, length(unique(names(grid.colB)))-1), 10*10, rep(0.05, length(unique(names(grid.colA)))-1), 10*10), start.degree = 220)
# now, the image with rotated labels
chordDiagram(df_chord_invitro, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col, 
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)
  #circos.axis(h = "top", labels.cex = 0.5, sector.index = sector.name, track.index = 2, major.tick = FALSE, labels = FALSE)
}, bg.border = NA)
circos.clear()
dev.off()

# 4. READ and REWRITE the SVG coordinate system to force the oval
svg_text <- readLines(paste0(FIGURE_PATH,"\\Fig3b_BipartiteGraph.svg"))
# Find where the main content group starts and apply a horizontal stretch factor (e.g., scale(1.6, 1))
# This explicitly alters the SVG matrix math to turn the circular paths into ellipses
svg_text <- gsub("<g id=", '<g transform="scale(1.0, 1.5)" id=', svg_text)
# Write the final modified oval file
writeLines(svg_text, paste0(FIGURE_PATH,"\\Fig3b_BipartiteGraph.svg"))













#ChordDiagram (invivo)
plot_data_invivo = plot_data[rowMeans(plot_data[,idx_invivo_efficacy_plot_data]) >= 5,]

final_label_with_newline = plot_data_invivo$EvenShorterLabel
to_label_modified = sapply(1:nrow(plot_data_invivo),function(i){
  if (plot_data_invivo$CompoundClass[i] == "Peptides"){
    return(paste0("Peptide: ", plot_data_invivo$Compound[i]))
  }
  return(gsub("-[(].*","", plot_data_invivo$Compound[i]))
})

df_chord_invivo = data.frame(from=to_label_modified, to=final_label_with_newline, value=rowMeans(plot_data_invivo[,idx_invivo_efficacy_plot_data]))

#grid.col <- setNames(rainbow(length(c(df_chord_invitro$from, df_chord_invitro$to))), unique(c(df_chord_invitro$from, df_chord_invitro$to)))
compound_class_color = c("#8DD3C7","#FFFFB3","#BEBADA","#FDB462","#DE9AEA")
names(compound_class_color) = c("Proteins","Lipids","Metabolites","Peptides","Oligosaccharides")
grid.colA = brewer.pal(6,"Set2")[as.numeric(factor(plot_data_invivo$FinalClass))]
names(grid.colA) = final_label_with_newline
grid.colB = compound_class_color[as.character(plot_data_invivo$CompoundClass)]
names(grid.colB) = to_label_modified
grid.col = c(grid.colA, grid.colB)




svg(paste0(FIGURE_PATH,"\\Fig3b_BipartiteGraph.svg"), width = 5, height = 5)
par(mar = c(3, 6, 4, 6), mfrow = c(1, 1), xpd = NA)
circos.par(gap.after = c(rep(2, length(unique(names(grid.colB)))-1), 10, rep(2, length(unique(names(grid.colA)))-1), 10))
# now, the image with rotated labels
chordDiagram(df_chord_invivo, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col,
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)
  #circos.axis(h = "top", labels.cex = 0.5, sector.index = sector.name, track.index = 2, major.tick = FALSE, labels = FALSE)
}, bg.border = NA)
dev.off()




