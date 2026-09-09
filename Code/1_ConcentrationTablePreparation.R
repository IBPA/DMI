#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\FilePath.R"))

#Molecule Data (6,718 molecules)
molecule_data = read.csv(paste0(DATA_PATH,"\\Molecules.csv"), header = T, check.names = F)

#Raw Conc Data (6,718 molecules)
raw_concentration_data = read.csv(paste0(DATA_PATH,"\\Concentration.csv"), header = T, check.names = F)

#Estimation Conc Data (4,434 molecules)
estimated_concentration_data = read.csv(paste0(DATA_PATH,"\\Estimated_Concentration.csv"), header = TRUE, check.names = F)

#Function for generating the concentration matrix
get_conc_matrix = function(conc_data, field = "Concentration Value"){
  conc_matrix = matrix(NA, nrow=length(unique(conc_data$`DMD Sample ID`)),
                       ncol=length(unique(conc_data$`DMD Molecule ID`)))
  rownames(conc_matrix) = unique(conc_data$`DMD Sample ID`)
  colnames(conc_matrix) = unique(conc_data$`DMD Molecule ID`)
  
  for (i in 1:nrow(conc_data)) {
    sample_id = conc_data$`DMD Sample ID`[i]
    molecule_id = conc_data$`DMD Molecule ID`[i]
    concentration_value = conc_data[i,field]
    conc_matrix[sample_id, molecule_id] = concentration_value
  }
  return(conc_matrix)
}

get_conc_matrix_avg = function(conc_data, field = "Concentration Value"){
  #Assume zero concentration if it is not detected
  avg_conc_matrix = matrix(0, nrow=length(unique(conc_data$`DMD Product ID`)),
                           ncol=length(unique(conc_data$`DMD Molecule ID`)))
  rownames(avg_conc_matrix) = unique(conc_data$`DMD Product ID`)
  colnames(avg_conc_matrix) = unique(conc_data$`DMD Molecule ID`)
  
  for (product_id in unique(conc_data$`DMD Product ID`)) {
    sample_ids = unique(conc_data$`DMD Sample ID`[conc_data$`DMD Product ID` == product_id])
    for (molecule_id in unique(conc_data$`DMD Molecule ID`)) {
      concentration_values = conc_data[
        conc_data$`DMD Product ID` == product_id &
          conc_data$`DMD Molecule ID` == molecule_id,
        field
      ]
      filled_concentration_values = concentration_values
      filled_concentration_values[is.na(filled_concentration_values)] = 0
      avg_conc_matrix[product_id, molecule_id] = mean(filled_concentration_values)
    }
  }
  return(avg_conc_matrix)
}

#1. Absolute Abundance Data (876 Compounds)
raw_concentration_data_absolute_abundance = raw_concentration_data[which(raw_concentration_data$`Concentration Type`=="Absolute Abundance"), ]
concentration_matrix_absolute_abundnace = get_conc_matrix(raw_concentration_data_absolute_abundance, field = "Concentration Value Alt")
concentration_matrix_absolute_abundnace_avg = get_conc_matrix_avg(raw_concentration_data_absolute_abundance, field = "Concentration Value Alt")

#2. Estimated Abundance Data (4,434 Compounds)
concentration_matrix_estimated_concentration = get_conc_matrix(estimated_concentration_data, field="Estimated Concentration Value")
concentration_matrix_estimated_concentration_avg = get_conc_matrix_avg(estimated_concentration_data, field="Estimated Concentration Value")

#3. Ion Count Data (980 Compounds)
raw_concentration_data_ion_count = raw_concentration_data[which(raw_concentration_data$`Concentration Type`=="Ion Count"), ]
raw_concentration_data_ion_count = raw_concentration_data_ion_count[which(!(raw_concentration_data_ion_count$`DMD Molecule ID` %in% raw_concentration_data_absolute_abundance$`DMD Molecule ID`)), ]
concentration_matrix_ion_count = get_conc_matrix(raw_concentration_data_ion_count)
concentration_matrix_ion_count_avg = get_conc_matrix_avg(raw_concentration_data_ion_count)

#4. miRNA Data (428 RNAs)
raw_concentration_data_miRNA = raw_concentration_data[which(raw_concentration_data$`Concentration Type`=="Molecule Count"), ]
concentration_matrix_miRNA = get_conc_matrix(raw_concentration_data_miRNA)
concentration_matrix_miRNA_avg = get_conc_matrix_avg(raw_concentration_data_miRNA)

concentration_matrix_absolute_abundnace_avg[is.na(concentration_matrix_absolute_abundnace_avg)] = 0
concentration_matrix_estimated_concentration_avg[is.na(concentration_matrix_estimated_concentration_avg)] = 0
concentration_matrix_ion_count_avg[is.na(concentration_matrix_ion_count_avg)] = 0
concentration_matrix_miRNA_avg[is.na(concentration_matrix_miRNA_avg)] = 0

#Save the concentration matrices

save(concentration_matrix_absolute_abundnace,
     concentration_matrix_estimated_concentration,
     concentration_matrix_ion_count,
     concentration_matrix_miRNA,
     concentration_matrix_absolute_abundnace_avg,
     concentration_matrix_estimated_concentration_avg,
     concentration_matrix_ion_count_avg,
     concentration_matrix_miRNA_avg,
     file = paste0(PROCESSED_DATA_PATH,"\\ConcentrationMatrices.rdata"))

data_pca_tsne = cbind(concentration_matrix_absolute_abundnace_avg, concentration_matrix_estimated_concentration_avg)
save(data_pca_tsne, file=paste0(PROCESSED_DATA_PATH,"\\compound_with_concentration.rdata"))