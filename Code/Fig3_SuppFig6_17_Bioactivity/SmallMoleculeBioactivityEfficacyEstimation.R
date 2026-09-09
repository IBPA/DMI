#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\BioassayDescription_Selected_Manual.csv")) || 
    !file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivity_FittedCurve.csv"))){
  stop("Please run the script 'GetSmallMoleculeBioactivity.r' first to generate the processed data files and prepare the manually selected bioassays (BioassayDescription_Selected_Manual.csv).")
}

#Select the target curves
aid_description_selected = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\BioassayDescription_Selected_Manual.csv"))
merged_aid_data_fitted_curve_addl_labeled = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivity_FittedCurve.csv"))
merged_aid_data_fitted_curve_addl_labeled_selected = merged_aid_data_fitted_curve_addl_labeled[merged_aid_data_fitted_curve_addl_labeled$AID %in% aid_description_selected$AID,]
merged_aid_data_fitted_curve_addl_labeled_selected = merged_aid_data_fitted_curve_addl_labeled_selected[merged_aid_data_fitted_curve_addl_labeled_selected$Efficacy != 0,]


get_efficacy = function(curve_data){
  conc_matrix = t(do.call(cbind, lapply(merged_aid_data_fitted_curve_addl_labeled_selected$DMD_ID, function(x){data_pca_tsne[,x,drop=F]})))
  colnames(conc_matrix) = rownames(data_pca_tsne)
  
  efficacy_matrix = sapply(1 : nrow(curve_data), function(i){
    hill_ac50 = (10^curve_data$Fit_LogAC50[i])*1e6 #Convert to μM
    hill_slope = curve_data$Fit_HillSlope[i]
    hill_min = curve_data$Fit_ZeroActivity[i]
    hill_max = curve_data$Fit_InfiniteActivity[i]
    
    mol_data_idx = which(mol_data$`DMD ID` == curve_data$DMD_ID[i])
    mw = mol_data$`Molecular Weight`[mol_data_idx]
    cur_conc_data = (conc_matrix[i,]/240*1000)/mw #ug/240mL to uM
    cur_conc_data_serum = cur_conc_data * (0.02/136)  #In vivo concentration adjustment
    
    cur_par = c(hill_ac50, hill_slope, hill_min, hill_max)
    
    efficacy_values = evalHillModel(cur_conc_data, cur_par)
    efficacy_values = (efficacy_values - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    
    efficacy_values_invivo = evalHillModel(cur_conc_data_serum, cur_par)
    efficacy_values_invivo = (efficacy_values_invivo - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    
    efficacy_values_invivo_100x = evalHillModel(cur_conc_data_serum*100, cur_par)
    efficacy_values_invivo_100x = (efficacy_values_invivo_100x - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    efficacy_values_invivo_10x = evalHillModel(cur_conc_data_serum*10, cur_par)
    efficacy_values_invivo_10x = (efficacy_values_invivo_10x - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    efficacy_values_invivo_01x = evalHillModel(cur_conc_data_serum*0.1, cur_par)
    efficacy_values_invivo_01x = (efficacy_values_invivo_01x - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    efficacy_values_invivo_001x = evalHillModel(cur_conc_data_serum*0.01, cur_par)
    efficacy_values_invivo_001x = (efficacy_values_invivo_001x - hill_min)/(hill_max - hill_min)*curve_data$Efficacy[i]
    
    mean_efficacy_value_invivo_100x = mean(efficacy_values_invivo_100x)
    mean_efficacy_value_invivo_10x = mean(efficacy_values_invivo_10x)
    mean_efficacy_value_invivo_01x = mean(efficacy_values_invivo_01x)
    mean_efficacy_value_invivo_001x = mean(efficacy_values_invivo_001x)
    
    return(c(cur_conc_data, cur_conc_data_serum, efficacy_values, efficacy_values_invivo,
             mean_efficacy_value_invivo_001x, mean_efficacy_value_invivo_01x, mean_efficacy_value_invivo_10x, mean_efficacy_value_invivo_100x))
  })
  efficacy_matrix = t(efficacy_matrix)
  colnames(efficacy_matrix) = c(rep(colnames(conc_matrix),4), paste0("Avg_Invivo-",c(0.01,0.1,10,100),"x"))
  
  return(efficacy_matrix)
}

efficacy_results = get_efficacy(merged_aid_data_fitted_curve_addl_labeled_selected)

write.csv(cbind(merged_aid_data_fitted_curve_addl_labeled_selected, efficacy_results, 
                rowMeans(efficacy_results[,121:180]), rowMeans(efficacy_results[,181:240]),
                sapply(merged_aid_data_fitted_curve_addl_labeled_selected$DMD_ID, map_compound_name)), 
          paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"), row.names = F) #The final milk small molecule bioactivity summary file, including the predicted in vivo efficacy values


#Note: Additional manual labeling is required for the next step (see MilkSmallMoleculeBioactivitySummary_Manual.csv as reference)