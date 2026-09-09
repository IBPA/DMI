#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

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
  Brand = product_data$`Brand Code`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)],
  Opacity = product_data$`Packaging Opaqueness`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)],
  Temperature = product_data$`Local Temperature`[match(rownames(ordered_compound_percentage), product_data$`DMD ID`)],
  DayToExpiration = difftime(as.POSIXct(product_data$`Expiration Date`,format="%Y/%m/%d"),as.POSIXct(product_data$`Purchase Date`,format="%Y/%m/%d %H:%M"))
)
rownames(annotation) = rownames(annotation)
annotation$Brand = factor(paste0("Brand ",annotation$Region, c("A","B")[annotation$Brand]), 
                          levels=paste0("Brand ",sapply(1:10,function(x){rep(x,2)}), rep(c("A","B"),10)))
annotation$Region = factor(regions[annotation$Region],levels=regions)
annotation$Time = factor(c("T1","T2","T3")[annotation$Time], levels=c("T1","T2","T3"))
annotation$Opacity = factor(annotation$Opacity, levels=c("Clear","Opaque"))


ordered_compound_percentage_nonzero = ordered_compound_percentage[,colSums(ordered_compound_percentage) > 0]
p_temp <- sapply(colnames(ordered_compound_percentage_nonzero), function(x){
  df = data.frame(
    Value = ordered_compound_percentage_nonzero[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Opacity, data=df)  # removes Temp + interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_temp <- p.adjust(p_temp, method="BH")

p_opacity <- sapply(colnames(ordered_compound_percentage_nonzero), function(x){
  df = data.frame(
    Value = ordered_compound_percentage_nonzero[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Temperature, data=df)  # removes Opacity + interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_opacity <- p.adjust(p_opacity, method="BH")

p_interaction <- sapply(colnames(ordered_compound_percentage_nonzero), function(x){
  df = data.frame(
    Value = ordered_compound_percentage_nonzero[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Opacity + Temperature, data=df)  # removes interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_interaction <- p.adjust(p_interaction, method="BH")


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

p_temp_molecular_class = lapply(colnames(sum_molecular_class), function(x){
  df = data.frame(
    Value = sum_molecular_class[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  df = df[which(!is.na(df$DayToExpiration)),]
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Opacity, data=df)  # removes Temp + interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_temp_molecular_class <- p.adjust(p_temp_molecular_class, method="BH")

p_opacity_molecular_class = lapply(colnames(sum_molecular_class), function(x){
  df = data.frame(
    Value = sum_molecular_class[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  df = df[which(!is.na(df$DayToExpiration)),]
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Temperature, data=df)  # removes Opacity + interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_opacity_molecular_class <- p.adjust(p_opacity_molecular_class, method="BH")

p_interaction_molecular_class = lapply(colnames(sum_molecular_class), function(x){
  df = data.frame(
    Value = sum_molecular_class[,x],
    Region = annotation$Region,
    Time = annotation$Time,
    Brand = annotation$Brand,
    Opacity = annotation$Opacity,
    Temperature = annotation$Temperature - mean(annotation$Temperature), #Centering the temperature for better model convergence
    DayToExpiration = annotation$DayToExpiration - mean(annotation$DayToExpiration, na.rm=T)
  )
  df = df[which(!is.na(df$DayToExpiration)),]
  
  fit_full <- lm(Value ~ Opacity * Temperature, data=df)
  fit_red  <- lm(Value ~ Opacity + Temperature, data=df)  # removes interaction
  anova(fit_red, fit_full)$`Pr(>F)`[2]
})
q_interaction_molecular_class <- p.adjust(p_interaction_molecular_class, method="BH")



cor_val_all = sapply(1:ncol(ordered_compound_percentage), function(i){
  cor(as.numeric(ordered_compound_percentage[,i]), as.numeric(annotation$Temperature), method = "pearson")
})
cor_pval_all = sapply(1:ncol(ordered_compound_percentage), function(i){
  cor.test(as.numeric(ordered_compound_percentage[,i]), as.numeric(annotation$Temperature), method = "pearson")$p.value
})
q_cor_pval_all = p.adjust(cor_pval_all, method = "BH")
spearman_cor_val_all = sapply(1:ncol(ordered_compound_percentage), function(i){
  cor(as.numeric(ordered_compound_percentage[,i]), as.numeric(annotation$Temperature), method = "spearman")
})
spearman_cor_pval_all = sapply(1:ncol(ordered_compound_percentage), function(i){
  cor.test(as.numeric(ordered_compound_percentage[,i]), as.numeric(annotation$Temperature), method = "spearman")$p.value
})
q_spearman_cor_pval_all = p.adjust(spearman_cor_pval_all, method = "BH")


cor_val_opaque = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")
})
cor_pval_opaque = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor.test(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")$p.value
})
q_cor_pval_opaque = p.adjust(cor_pval_opaque, method = "BH")
spearman_cor_val_opaque = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")
})
spearman_cor_pval_opaque = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor.test(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")$p.value
})
q_spearman_cor_pval_opaque = p.adjust(spearman_cor_pval_opaque, method = "BH")

cor_val_clear = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")
})
cor_pval_clear = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor.test(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")$p.value
})
q_cor_pval_clear = p.adjust(cor_pval_clear, method = "BH")
spearman_cor_val_clear = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")
})
spearman_cor_pval_clear = sapply(1:ncol(ordered_compound_percentage), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor.test(as.numeric(ordered_compound_percentage[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")$p.value
})
q_spearman_cor_pval_clear = p.adjust(spearman_cor_pval_clear, method = "BH")



cor_val_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  cor(as.numeric(sum_molecular_class[,i]), as.numeric(annotation$Temperature), method = "pearson")
})
cor_pval_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  cor.test(as.numeric(sum_molecular_class[,i]), as.numeric(annotation$Temperature), method = "pearson")$p.value
})
q_cor_pval_molecule_class = p.adjust(cor_pval_molecule_class, method = "BH")
spearman_cor_val_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  cor(as.numeric(sum_molecular_class[,i]), as.numeric(annotation$Temperature), method = "spearman")
})
spearman_cor_pval_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  cor.test(as.numeric(sum_molecular_class[,i]), as.numeric(annotation$Temperature), method = "spearman")$p.value
})
q_spearman_cor_pval_molecule_class = p.adjust(spearman_cor_pval_molecule_class, method = "BH")

cor_val_opaque_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")
})
cor_pval_opaque_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor.test(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")$p.value
})
q_cor_pval_opaque_molecule_class = p.adjust(cor_pval_opaque_molecule_class, method = "BH")
spearman_cor_val_opaque_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")
})
spearman_cor_pval_opaque_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Opaque")
  cor.test(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")$p.value
})
q_spearman_cor_pval_opaque_molecule_class = p.adjust(spearman_cor_pval_opaque_molecule_class, method = "BH")

cor_val_clear_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")
})
cor_pval_clear_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor.test(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "pearson")$p.value
})
q_cor_pval_clear_molecule_class = p.adjust(cor_pval_clear_molecule_class, method = "BH")
spearman_cor_val_clear_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")
})
spearman_cor_pval_clear_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  idx = which(annotation$Opacity == "Clear")
  cor.test(as.numeric(sum_molecular_class[idx,i]), as.numeric(annotation$Temperature[idx]), method = "spearman")$p.value
})
q_spearman_cor_pval_clear_molecule_class = p.adjust(spearman_cor_pval_clear_molecule_class, method = "BH")


pval_ttest_opacity = sapply(1:ncol(ordered_compound_percentage), function(i){
  t.test(ordered_compound_percentage[which(annotation$Opacity == "Clear"), i], ordered_compound_percentage[which(annotation$Opacity == "Opaque"), i])$p.value
})
qval_ttest_opacity = p.adjust(pval_ttest_opacity, method = "BH")

pval_ttest_opacity_molecule_class = sapply(1:ncol(sum_molecular_class), function(i){
  t.test(sum_molecular_class[which(annotation$Opacity == "Clear"), i], sum_molecular_class[which(annotation$Opacity == "Opaque"), i])$p.value
})
qval_ttest_opacity_molecule_class = p.adjust(pval_ttest_opacity_molecule_class, method = "BH")



#Final significant items
significant_compound_temp_idx = which(q_temp < 0.1 & q_cor_pval_all[1:5172] < 0.1 & q_spearman_cor_pval_all[1:5172] < 0.1)
significant_compound_opacity_idx = which(q_opacity < 0.1 & qval_ttest_opacity[1:5172] < 0.1)
#significant_interaction_idx = which(q_interaction < 0.1 & 
#                                      ((q_cor_pval_clear[1:5172] < 0.1 & q_spearman_cor_pval_clear[1:5172] < 0.1) | 
#                                      (q_cor_pval_opaque[1:5172] < 0.1 & q_spearman_cor_pval_opaque[1:5172] < 0.1)))

#significant_interaction_temp_class_idx = which(q_temp_molecular_class < 0.1 & q_cor_pval_molecule_class < 0.1 & q_spearman_cor_pval_molecule_class < 0.1)
#significant_interaction_opacity_class_idx = which(q_opacity_molecular_class < 0.1 & qval_ttest_opacity_molecule_class < 0.1)
significant_interaction_class_idx = which(q_interaction_molecular_class < 0.1 & 
                                      ((q_cor_pval_opaque_molecule_class < 0.1 & q_spearman_cor_pval_opaque_molecule_class < 0.1) | 
                                      (q_cor_pval_clear_molecule_class < 0.1 & q_spearman_cor_pval_clear_molecule_class < 0.1)))




#Pos temperature correlation
df_DEELNKLLG = data.frame(
  Percentage = ordered_compound_percentage[, significant_compound_temp_idx[4]]*100,
  Temperature = sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Local Temperature"]}),
  Region = factor(regions[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)

svg(paste0(FIGURE_PATH,"\\SuppFig9a_NegCorr_Temp.svg"), width = 3.5, height = 3.5)
p_corr = ggplot(df_DEELNKLLG, aes(y=Percentage, x=Temperature, alpha=Time, shape = Brand, fill = Region)) +
  geom_point(size=3, stroke=1.0) +
  theme_classic() +
  labs(title="", x="Purchase location temperature (°F)", y="Content (Peptide:DEELNKLLG) (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
  guides(fill=guide_legend(override.aes=list(shape=21)), 
         alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
         shape=guide_legend(override.aes=list(alpha=1)))
print(p_corr)
dev.off()

#Neg temperature correlation
df_LHQKHQRAKRAVSHEDQFLRLDF = data.frame(
  Percentage = ordered_compound_percentage[, significant_compound_temp_idx[3]]*100,
  Temperature = sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Local Temperature"]}),
  Region = factor(regions[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)

svg(paste0(FIGURE_PATH,"\\SuppFig9b_PosCorr_Temp.svg"), width = 4.5, height = 4.5)
p_corr = ggplot(df_LHQKHQRAKRAVSHEDQFLRLDF, aes(y=Percentage, x=Temperature, alpha=Time, shape = Brand, fill = Region)) +
  geom_point(size=3, stroke=1.0) +
  theme_classic() +
  labs(title="", x="Purchase location temperature (°F)", y="Content (Peptide:LHQKHQRAKRAVSHEDQFLRLDF) (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "right") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
  guides(fill=guide_legend(override.aes=list(shape=21)), 
         alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
         shape=guide_legend(override.aes=list(alpha=1)))
print(p_corr)
dev.off()


#Boxplot (opaque)
df_packaging_opaqueness = data.frame(
  Percentage = unlist(lapply(significant_compound_opacity_idx[1], function(mol_id){
    return(ordered_compound_percentage[, mol_id]*100)
  })),
  Compound = factor(rep(sapply(names(significant_compound_opacity_idx)[1], map_compound_name), each = nrow(ordered_compound_percentage)),
                    levels = sapply(names(significant_compound_opacity_idx)[1], map_compound_name)),
  PackagingOpaqueness = factor(rep(product_data[sapply(rownames(ordered_compound_percentage),function(x){which(x==product_data$`DMD ID`)}),"Packaging Opaqueness"],
                                   times = length(significant_compound_opacity_idx[1])),
                               levels = c("Clear","Opaque"))
)
svg(paste0(FIGURE_PATH,"\\SuppFig9c_Opaqueness.svg"), width = 2, height = 3.5)
p_packaging_all = ggplot(df_packaging_opaqueness, aes(y=Percentage, x=Compound, fill=PackagingOpaqueness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.5, size=1.0, position=position_jitterdodge(jitter.width=0.25, dodge.width=0.75)) +
  theme_classic() +
  labs(title="", x="", y="Content (Aspartic acid) (% w/w)",fill="Packaging\nOpaqueness") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_fill_manual(values = brewer.pal(2,"Set2"))
print(p_packaging_all)
dev.off()


#Correlation (Lipid)
df_lipid = data.frame(
  Percentage = sum_molecular_class[, significant_interaction_class_idx[1]]*100,
  Temperature = sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Local Temperature"]}),
  Region = factor(regions[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Packaging = factor(sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Packaging Opaqueness"]}), levels = c("Clear","Opaque")),
  Brand = c("Brand A","Brand B")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(ordered_compound_percentage),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)

svg(paste0(FIGURE_PATH,"\\SuppFig9d_Corr_Clear.svg"), width = 3.5, height = 3.5)
p_corr = ggplot(df_lipid[which(df_lipid$Packaging=="Clear"),], aes(y=Percentage, x=Temperature, alpha=Time, shape = Brand, fill = Region)) +
  geom_point(size=3, stroke=1.0) +
  theme_classic() +
  labs(title="", x="Purchase location temperature (°F)", y="Content (Total Lipid) (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")[c(1,2,3,4,5,7,8,9)]) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
  guides(fill=guide_legend(override.aes=list(shape=21)), 
         alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
         shape=guide_legend(override.aes=list(alpha=1)))
print(p_corr)
dev.off()

svg(paste0(FIGURE_PATH,"\\SuppFig9e_Corr_Opaque.svg"), width = 3.5, height = 3.5)
p_corr = ggplot(df_lipid[which(df_lipid$Packaging=="Opaque"),], aes(y=Percentage, x=Temperature, alpha=Time, shape = Brand, fill = Region)) +
  geom_point(size=3, stroke=1.0) +
  theme_classic() +
  labs(title="", x="Purchase location temperature (°F)", y="Content (Total Lipid) (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")[c(2,5,6,8,10)]) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
  guides(fill=guide_legend(override.aes=list(shape=21)), 
         alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
         shape=guide_legend(override.aes=list(alpha=1)))
print(p_corr)
dev.off()

