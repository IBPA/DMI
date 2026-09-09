#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

if (!file.exists(paste0(DATA_PATH,"\\Bioactivity\\ManuallyCuratedPeptideIC50Data.csv"))){ 
  stop("A manually curated peptide with IC50 values is needed. (ManuallyCuratedPeptideIC50Data.csv).")
}
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"))){ 
    stop("A representative Hill curve evaluated based on the small molecules dataset is needed. Please run the script 'SmallMoleculeBioactivityEfficacyEstimation.r' first to extract the Hill curves (MilkSmallMoleculeBioactivitySummary.csv).")
}


library(ggrastr)

curve_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"))
curve_data = curve_data[which(curve_data[,"Fit_ZeroActivity"] != curve_data[,"Fit_InfiniteActivity"]),]
fold_change = 10^seq(-5,5,0.1)
efficacy_curves = sapply(1 : nrow(curve_data), function(i){
  hpar = as.numeric(curve_data[i,c("Fit_LogAC50","Fit_HillSlope","Fit_ZeroActivity","Fit_InfiniteActivity")])
  hpar[1] = 10^(hpar[1])*1e6 #Convert logAC50 to AC50
  
  efficacy_values = basicdrm::evalHillModel(fold_change*(hpar[1]), hpar)
  efficacy_values = (efficacy_values - hpar[3])/(hpar[4] - hpar[3])*100
  
  return(efficacy_values)
})
efficacy_curves = t(efficacy_curves)
ratio_efficacy_eval_mean = apply(efficacy_curves,2,mean)
ratio_efficacy_eval_std = apply(efficacy_curves,2,sd)

ratio_efficacy_eval_mean_curve = fitHillModel(
  10^seq(-5,5,0.1),
  as.numeric(ratio_efficacy_eval_mean),
  c(1,2),
  start=c(1,2,0,100),
)
#plot(evalHillModel(10^seq(-5,5,0.1), ratio_efficacy_eval_mean_curve$coefficients) - ratio_efficacy_eval_mean) #Max fitted error: 3.17(%)



svg(paste0(FIGURE_PATH,"\\SuppFig17a_HillCurves.svg"), width = 4, height = 4)
df_efficacy_plot = data.frame(
  FoldChange = rep(seq(-5,5,0.1), each = nrow(curve_data)),
  Efficacy = as.vector(efficacy_curves)
)
df_curve = data.frame(
  x = seq(-5,5,0.1),
  y_mean = (ratio_efficacy_eval_mean),
  y_upper = (ratio_efficacy_eval_mean + 2*ratio_efficacy_eval_std),
  y_lower = (ratio_efficacy_eval_mean - 2*ratio_efficacy_eval_std)
)
h = ggplot(df_efficacy_plot) +
  geom_point_rast(aes(x=FoldChange, y=Efficacy), alpha=0.1, color="grey") +
  theme_classic() +
  labs(title="", x="Log10(Concentration/IC50)", y="Normalized Efficacy (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(-5,5)) +
  ylim(c(-20,120)) + 
  geom_ribbon(data = df_curve, aes(x = x, ymin = y_lower, ymax = y_upper), fill="#AAEEAA", alpha=0.3) + 
  geom_line(data=df_curve, aes(x=x, y=y_mean),color="black")
print(h)
dev.off()

protein_peptide_data = read.csv(paste0(DATA_PATH,"\\Bioactivity\\ManuallyCuratedPeptideIC50Data.csv"), header = T, check.names = F)
protein_peptide_data$DMD_ID_Peptide = sapply(protein_peptide_data$Peptide, function(x){
  idx = which(x == mol_data$`Chemical Composition`)
  if (length(idx) == 0) return(NA)
  if (length(idx) == 1) return(mol_data$`DMD ID`[idx])
  print(x)
  stop('!')
})
protein_peptide_data$n_intervals = sapply(strsplit(protein_peptide_data$Intervals,","),length)

get_efficacy = function(protein_peptide_data){
  
  efficacy_matrix = sapply(1 : nrow(protein_peptide_data), function(i){
    if (is.na(protein_peptide_data$DMD_ID_Peptide[i])){
      cur_conc_data_peptide = rep(0, nrow(data_pca_tsne))
    }else{
      mol_data_idx = which(mol_data$`DMD ID` == protein_peptide_data$DMD_ID_Peptide[i])
      mw = mol_data$`Molecular Weight`[mol_data_idx]
      cur_conc_data_peptide = (data_pca_tsne[,protein_peptide_data$DMD_ID_Peptide[i]]/240*1000)/mw #ug/240mL to uM
    }
    mol_data_idx = which(mol_data$`DMD ID` == protein_peptide_data$DMD_ID_Protein[i])
    mw = mol_data$`Molecular Weight`[mol_data_idx]
    cur_conc_data_protein = (data_pca_tsne[,protein_peptide_data$DMD_ID_Protein[i]]/240*1000)/mw #ug/240mL to uM
    
    cur_conc_data = cur_conc_data_protein * protein_peptide_data$n_intervals[i] + cur_conc_data_peptide
    cur_conc_data_serum = cur_conc_data * 0.0006/100
    efficacy_values = evalHillModel(cur_conc_data/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients)
    efficacy_values_invivo = evalHillModel(cur_conc_data_serum/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients)
    
    mean_efficacy_value_invivo_100x = mean(evalHillModel(100*cur_conc_data_serum/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients))
    mean_efficacy_value_invivo_10x = mean(evalHillModel(10*cur_conc_data_serum/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients))
    mean_efficacy_value_invivo_01x = mean(evalHillModel(0.1*cur_conc_data_serum/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients))
    mean_efficacy_value_invivo_001x = mean(evalHillModel(0.01*cur_conc_data_serum/protein_peptide_data$`Final IC50 (uM)`[i], ratio_efficacy_eval_mean_curve$coefficients))
    
    return(c(cur_conc_data, cur_conc_data_serum, efficacy_values, efficacy_values_invivo,
             cur_conc_data_peptide, cur_conc_data_protein, 
             mean_efficacy_value_invivo_001x, mean_efficacy_value_invivo_01x, mean_efficacy_value_invivo_10x, mean_efficacy_value_invivo_100x
             ))
  })
  efficacy_matrix = t(efficacy_matrix)
  colnames(efficacy_matrix) = c(rep(rownames(data_pca_tsne),6), paste0("Avg_Invivo-",c(0.01,0.1,10,100),"x"))
  
  return(efficacy_matrix)
}

efficacy_results = get_efficacy(protein_peptide_data)

write.csv(cbind(protein_peptide_data, efficacy_results, 
                rowMeans(efficacy_results[,121:180]), rowMeans(efficacy_results[,181:240]),
                rowMeans(efficacy_results[,241:300]), rowMeans(efficacy_results[,301:360])), 
          paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkPeptideBioactivitySummary.csv"), row.names = F)


