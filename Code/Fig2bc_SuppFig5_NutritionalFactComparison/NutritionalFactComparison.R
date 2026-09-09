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

#Nutritional Facts Extraction:
nutritional_fact = sapply(strsplit(product_data$`Nutritional Facts`,","), function(x){
  result_val = sapply(x, function(y){
    y = gsub("\"","",y)
    y = gsub("{","",y,fixed=T)
    y = gsub("}","",y,fixed=T)
    kv = strsplit(y,":")[[1]]
    value = as.numeric(kv[2])
    return(value)
  })
  new_name = gsub("\"","",x)
  new_name = gsub("{","",new_name,fixed=T)
  new_name = gsub("}","",new_name,fixed=T)
  names(result_val) = sapply(strsplit(new_name,":"),function(x){x[1]})
  return(result_val)
})

nutritional_fact_matrix = matrix(NA, ncol=length(unique(names(unlist(nutritional_fact)))), nrow=length(product_data$`DMD ID`))
rownames(nutritional_fact_matrix) = product_data$`DMD ID`
colnames(nutritional_fact_matrix) = unique(names(unlist(nutritional_fact)))

for (i in 1 : length(nutritional_fact)){
  nutritional_fact_matrix[i, names(nutritional_fact[[i]])] = nutritional_fact[[i]]
}

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
colnames(sum_molecular_class)[ncol(sum_molecular_class)] = "Saccharide"


#Diff. Boxplot (All)
df_diff_boxplot = data.frame(
  quantified_value = c(sum_molecular_class[,"Proteins"]*240*1.03, 
                       sum_molecular_class[,"Lipids"]*240*1.03, 
                       sum_molecular_class[,"Saccharide"]*240*1.03, 
                       data_pca_tsne[,"DMD306416"]*1.03/1e6*1000),
  nutritional_fact_value = c(nutritional_fact_matrix[rownames(sum_molecular_class),"Protein"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
                                 nutritional_fact_matrix[rownames(sum_molecular_class),"Total Fat"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
                                 nutritional_fact_matrix[rownames(sum_molecular_class),"Total Sugars"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
                                 nutritional_fact_matrix[rownames(sum_molecular_class),"Cholesterol"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240),

  Class = factor(c(rep("Total\nProtein", nrow(sum_molecular_class)),
                   rep("Total\nFat", nrow(sum_molecular_class)),
                   rep("Total\nSugar", nrow(sum_molecular_class)),
                   rep("Cholesterol", nrow(sum_molecular_class))),
                 levels = c("Total\nSugar","Total\nProtein","Total\nFat","Cholesterol"))
  
)

df_diff_boxplot$RelDiff = (df_diff_boxplot$quantified_value - df_diff_boxplot$nutritional_fact_value) / df_diff_boxplot$nutritional_fact_value * 100

#Diff distribution
df_stacked_bar = data.frame(
  mean_value = c(mean(sum_molecular_class[,"Proteins"]), 
                 mean(sum_molecular_class[,"Lipids"]), 
                 mean(sum_molecular_class[,"Saccharide"]),
                 mean(nutritional_fact_matrix[rownames(sum_molecular_class),"Protein"])/(240*1.03),
                 mean(nutritional_fact_matrix[rownames(sum_molecular_class),"Total Fat"]/(240*1.03)),
                 mean(nutritional_fact_matrix[rownames(sum_molecular_class),"Total Sugars"]/(240*1.03))),
  
  Class = factor(rep(c("Total Protein","Total Fat","Total Sugar"),2), levels=c("Total Protein","Total Fat","Total Sugar")),
  Type = factor(c(rep("Quantified",3), rep("Nutritional Fact",3)), levels=c("Quantified","Nutritional Fact"))
)

library(patchwork)
svg(paste0(FIGURE_PATH,"\\Fig2bc.svg"), width = 7, height = 3.5)
p_diff_boxplot = ggplot(df_diff_boxplot, aes(x=Class, y=RelDiff, fill=Class)) +
  geom_boxplot() +
  theme_classic() +
  labs(title="", x="", y="Relative Difference (%)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), axis.title=element_blank(), legend.position = "none") +
  scale_fill_manual(values = brewer.pal(n = 4, name = "Set2")[c(1,2,3,4)]) + 
  geom_hline(yintercept = 0, linetype="dashed", color="red") + 
  scale_y_continuous(limits = c(-90, 90), breaks = seq(-75,75,25)) + 
  geom_jitter(position=position_jitter(width=0.2), size=1, shape=21, fill="black", alpha=0.3)
#print(p_diff_boxplot)
#dev.off()

#svg("Figures\\Q4_NutritionalFact_Diff_StackedBar.svg", width = 3, height = 3)
p_stacked_bar = ggplot(df_stacked_bar, aes(x=Type, y=mean_value*100, fill=Class)) +
  geom_bar(stat="identity", position="stack") +
  theme_classic() +
  labs(title="", x="", y="Content (% w/w)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), axis.title.x=element_blank(), legend.position = "none") +
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3","#CEAAE2")) + 
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,2), labels=c("0%","2%","4%","6%","8%","10%","12%"), expand=c(0,0))
#print(p_stacked_bar)

p = p_diff_boxplot + p_stacked_bar + plot_layout(widths = c(2,1))
print(p)
dev.off()



#Scatter Plot (Fat)
df_fat = data.frame(
  FatPercentage = sum_molecular_class[,"Lipids"]*240*1.03,
  FatNutritionalFact = nutritional_fact_matrix[rownames(sum_molecular_class),"Total Fat"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
  Region = factor(regions[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)
#svg("Figures\\Q4_NutritionalFact_Correlation_Fat.svg", width = 2, height = 2)
p_fat = ggplot(df_fat, aes(y=FatPercentage, x=FatNutritionalFact, fill = Region)) +
  geom_point(aes(color=Region),size=3, shape=21, color="black") +
  theme_classic() +
  labs(title="", x="Fat (Nutritional Fact) (g/240 mL)", y="Total Quantified Fat (g/240 mL)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"),
        axis.title=element_blank(), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) +
  xlim(c(0,10)) + ylim(c(0,10)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red")
#print(p_fat)
#dev.off()

#Scatter Plot (Sugar)
df_sugar = data.frame(
  SugarPercentage = sum_molecular_class[,"Saccharide"]*240*1.03,
  SugarNutritionalFact = nutritional_fact_matrix[rownames(sum_molecular_class),"Total Sugars"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
  Region = factor(regions[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)
#svg("Figures\\Q4_NutritionalFact_Correlation_Sugar.svg", width = 2, height = 2)
p_sugar = ggplot(df_sugar, aes(y=SugarPercentage, x=SugarNutritionalFact, fill = Region)) +
  geom_point(aes(color=Region),size=3, shape=21, color="black") +
  theme_classic() +
  labs(title="", x="Sugar (Nutritional Fact) (g/240 mL)", y="Total Quantified Saccharide (Sugar) (g/240 mL)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"),
        axis.title=element_blank(), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) +
  xlim(c(8,15)) + ylim(c(8,15)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red")
#print(p_sugar)
#dev.off()

#Scatter Plot (Protein)
df_protein = data.frame(
  ProteinPercentage = sum_molecular_class[,"Proteins"]*240*1.03,
  ProteinNutritionalFact = nutritional_fact_matrix[rownames(sum_molecular_class),"Protein"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
  Region = factor(regions[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)
#svg("Figures\\Q4_NutritionalFact_Correlation_Protein.svg", width = 2, height = 2)
p_protein = ggplot(df_protein, aes(y=ProteinPercentage, x=ProteinNutritionalFact, fill = Region)) +
  geom_point(aes(color=Region),size=3, shape=21, color="black") +
  theme_classic() +
  labs(title="", x="Protein (Nutritional Fact) (g/240 mL)", y="Total Quantified Protein (g/240 mL)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"),
        axis.title=element_blank(), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) +
  xlim(c(7,15)) + ylim(c(7,15)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red")
#print(p_protein)
#dev.off()

#Scatter Plot (Cholesterol)
df_cholesterol = data.frame(
  CholesterolPercentage = data_pca_tsne[,"DMD306416"]*1.03/1e6*1000,
  CholesterolNutritionalFact = nutritional_fact_matrix[rownames(sum_molecular_class),"Cholesterol"]/(nutritional_fact_matrix[rownames(sum_molecular_class),"Serving Size"])*240,
  Region = factor(regions[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Region Code"]})],levels = regions),
  Brand = c("Brand A","Brand B")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Brand Code"]})],
  Time = c("T1","T2","T3")[sapply(rownames(sum_molecular_class),function(x){product_data[which(x==product_data$`DMD ID`),"Time Point Code"]})]
)
#svg("Figures\\Q4_NutritionalFact_Correlation_Cholesterol.svg", width = 2, height = 2)
p_cholesterol = ggplot(df_cholesterol, aes(y=CholesterolPercentage, x=CholesterolNutritionalFact, fill = Region)) +
  geom_point(aes(color=Region),size=3, shape=21, color="black") +
  theme_classic() +
  labs(title="", x="Cholesterol (Nutritional Fact) (mg/240 mL)", y="Total Quantified Cholesterol (mg/240 mL)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), 
        axis.title=element_blank(), legend.position = "none") + 
  scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) +
  xlim(c(0,35)) + ylim(c(0,35)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red")
#print(p_cholesterol)
#dev.off()

svg(paste0(FIGURE_PATH,"\\SuppFig5.svg"), width = 7, height = 2)
p_all = p_sugar + p_protein + p_fat + p_cholesterol + plot_layout(ncol=4)
print(p_all)
dev.off()
