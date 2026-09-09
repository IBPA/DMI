#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))


if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary_Manual.csv"))){
  stop("Please run the script 'SmallMoleculeBioactivityEfficacyEstimation.r' with necessary manual labeling first to estimate the milk small molecule bioactivity efficacy.")
}
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkPeptideBioactivitySummary.csv"))){
  stop("Please run the script 'PeptideEfficacyEstimation.r' first to estimate the milk peptide bioactivity efficacy.")
}


efficacy_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary_Manual.csv"))
n_column_metadata_efficacy = 27
idx_conc_data = n_column_metadata_efficacy + (1:60)
idx_invivo_conc_data = n_column_metadata_efficacy + (61:120)
idx_efficacy_data = n_column_metadata_efficacy + (121:180)
idx_invivo_efficacy_data = n_column_metadata_efficacy + (181:240)

peptide_efficacy_data = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkPeptideBioactivitySummary.csv"))
n_column_metadata_peptide_efficacy = 15
idx_conc_peptide_data = n_column_metadata_peptide_efficacy + (301:360)
idx_invivo_conc_peptide_data = n_column_metadata_peptide_efficacy + (61:120)
idx_peptide_efficacy_data = n_column_metadata_peptide_efficacy + (121:180)
idx_invivo_peptide_efficacy_data = n_column_metadata_peptide_efficacy + (181:240)


plot_efficacy_matrix = NULL
plot_efficacy_matrix_invivo = NULL
bioactivity_name = c()
corresponding_compound_name = c()
corresponding_compound_class = c()
average_compound_conc = c()
reference_info = c()
efficacy_invivo_001x = c()
efficacy_invivo_01x = c()
efficacy_invivo_10x = c()
efficacy_invivo_100x = c()

efficacy_data = efficacy_data[order(rowMeans(efficacy_data[,idx_invivo_efficacy_data]),decreasing=T),]
peptide_efficacy_data = peptide_efficacy_data[order(rowMeans(peptide_efficacy_data[,idx_invivo_peptide_efficacy_data]),decreasing=T),]

for (i in 1 : nrow(efficacy_data)){
  cur_bioactivity_name = efficacy_data$Purpose[i]
  if (!(cur_bioactivity_name %in% bioactivity_name)){
    if (mean(as.numeric(efficacy_data[i,idx_conc_data])) == 0) next
    
    bioactivity_name = c(bioactivity_name, cur_bioactivity_name)
    corresponding_compound_name = c(corresponding_compound_name, efficacy_data$CompoundName[i])
    average_compound_conc = c(average_compound_conc, mean(as.numeric(efficacy_data[i,idx_conc_data])))
    corresponding_compound_class = c(corresponding_compound_class, molecular_classification(efficacy_data$DMD_ID[i]))
    reference_info = c(reference_info, efficacy_data$AID[i])
    efficacy_invivo_001x = c(efficacy_invivo_001x, efficacy_data$Avg_Invivo.0.01x[i])
    efficacy_invivo_01x = c(efficacy_invivo_01x, efficacy_data$Avg_Invivo.0.1x[i])
    efficacy_invivo_10x = c(efficacy_invivo_10x, efficacy_data$Avg_Invivo.10x[i])
    efficacy_invivo_100x = c(efficacy_invivo_100x, efficacy_data$Avg_Invivo.100x[i])
    
    plot_efficacy_matrix = rbind(plot_efficacy_matrix, efficacy_data[i,idx_efficacy_data])
    plot_efficacy_matrix_invivo = rbind(plot_efficacy_matrix_invivo, efficacy_data[i,idx_invivo_efficacy_data])
  }
}
for (i in 1 : nrow(peptide_efficacy_data)){
  cur_bioactivity_name = peptide_efficacy_data$Function[i]
  if (!(cur_bioactivity_name %in% bioactivity_name)){
    if (mean(as.numeric(peptide_efficacy_data[i,idx_conc_peptide_data])) == 0) next
    
    bioactivity_name = c(bioactivity_name, cur_bioactivity_name)
    corresponding_compound_name = c(corresponding_compound_name, peptide_efficacy_data$Peptide[i])
    average_compound_conc = c(average_compound_conc, mean(as.numeric(peptide_efficacy_data[i,idx_conc_peptide_data])))
    corresponding_compound_class = c(corresponding_compound_class, "Peptides")
    reference_info = c(reference_info, peptide_efficacy_data$Selected.Literature[i])
    efficacy_invivo_001x = c(efficacy_invivo_001x, peptide_efficacy_data$Avg_Invivo.0.01x[i])
    efficacy_invivo_01x = c(efficacy_invivo_01x, peptide_efficacy_data$Avg_Invivo.0.1x[i])
    efficacy_invivo_10x = c(efficacy_invivo_10x, peptide_efficacy_data$Avg_Invivo.10x[i])
    efficacy_invivo_100x = c(efficacy_invivo_100x, peptide_efficacy_data$Avg_Invivo.100x[i])
    
    plot_efficacy_matrix = rbind(plot_efficacy_matrix, peptide_efficacy_data[i,idx_peptide_efficacy_data])
    plot_efficacy_matrix_invivo = rbind(plot_efficacy_matrix_invivo, peptide_efficacy_data[i,idx_invivo_peptide_efficacy_data])
  }
}

plot_matrix_summary = cbind(bioactivity_name, reference_info, corresponding_compound_name, average_compound_conc,
                            plot_efficacy_matrix, plot_efficacy_matrix_invivo,
                            efficacy_invivo_001x, efficacy_invivo_01x, efficacy_invivo_10x, efficacy_invivo_100x, corresponding_compound_class
                            )
colnames(plot_matrix_summary) = c("Bioactivity","AID/Study","Compound","Avg. Conc.",
                                  colnames(efficacy_data)[idx_conc_data], colnames(peptide_efficacy_data)[idx_conc_peptide_data],
                                  paste0("Avg_Invivo_",c(0.01,0.1,10,100),"x"), "CompoundClass")
write.csv(plot_matrix_summary, paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\Plot_Efficacy_Matrix_Summary.csv"))

#Additional name abbreviation and classification is required. Please check the manually modified file "Plot_Efficacy_Matrix_Summary_ManualLabel.csv" as reference.
