if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\FilePath.R"))

#Product Data
product_data = read.csv(paste0(DATA_PATH,"\\Products.csv"), header = T, check.names = F)
if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\compound_with_concentration.rdata"))){
  stop("Please run the script '1_ConcentrationTablePreparation.R' to generate the processed data first.")
}
load(paste0(PROCESSED_DATA_PATH,"\\compound_with_concentration.rdata"))
regions = c("Northeast","Mid-South","West","South Central","Plains","California-1","Southeast","Great Lakes","California-2","California-3")

molecular_classification = function(mol_id){
  if (is.na(mol_id)){
    return(NA)
  }
  mol_class = mol_data$`Molecule Classification`[which(mol_data$`DMD ID` == mol_id)]
  mol_lab = mol_data$`Omic Lab`[which(mol_data$`DMD ID` == mol_id)]
  
  if (grepl("Peptidomics", mol_lab)){
    return("Peptides")
  }else if (grepl("Proteomics", mol_lab)){
    return("Proteins")
  }else if (grepl("Lipidomics", mol_lab)){
    return("Lipids")
  }else if (grepl("Glycomics", mol_lab)){
    return("Oligosaccharides")
  }else if (grepl("Targeted Metabolomics", mol_lab)){
    return("Metabolites")
  }else if (grepl("Untargeted Metabolomics", mol_lab)){
    return("Metabolites (Untargeted)")
  }else{
    return("miRNA") #Should not happen except miRNA
  }
}

mol_data = read.csv(paste0(DATA_PATH,"\\Molecules.csv"), header = T, check.names = F)
map_compound_name = function(DMD_ID){
  name = mol_data$`Molecule Name`[which(mol_data$`DMD ID` == DMD_ID)]
  if(length(name) == 0){
    return(DMD_ID)
  } else {
    return(name)
  }
}

library(jsonlite)
map_compound_name2 = function(DMD_ID){
  idx = which(mol_data$`DMD ID` == DMD_ID)
  name = mol_data$`Molecule Name`[idx]
  
  if (grepl("Peptide",mol_data$`Molecule Classification`[idx])){
    name = paste0("Peptide:", mol_data$`Chemical Composition`[idx])
  }else if (grepl("Protein",mol_data$`Molecule Classification`[idx])){
    id_info = fromJSON(mol_data$`External Database IDs`[idx])
    name = paste0(name,"\n(Uniprot:", id_info$UniProt[1],")")
  }
  return(name)
}

all_compound_percentage = data_pca_tsne/1e6/240/1.03 #Convert to percentage - assuming density of 1.03 g/mL
#Percentage of top compounds
mean_all_compound_percentage = colMeans(all_compound_percentage)
mean_all_compound_percentage = mean_all_compound_percentage[order(mean_all_compound_percentage, decreasing = T)]

ordered_compound_percentage = data_pca_tsne[, names(mean_all_compound_percentage)]/1e6/240/1.03 #Convert to percentage - assuming density of 1.03 g/mL
sum_top_N_compound_percentage = sapply(c(10,20,50,100,ncol(ordered_compound_percentage)), function(N){
  return(rowSums(ordered_compound_percentage[,1:N]))
})
colnames(sum_top_N_compound_percentage) = c("top_10", "top_20", "top_50", "top_100", "all_compounds")

#Wholemilk sample only
ordered_compound_percentage_wholemilk = ordered_compound_percentage[which(!rownames(ordered_compound_percentage) %in% c("DMD100028","DMD100030")),]
sum_top_N_compound_percentage_wholemilk = sapply(c(10,20,50,100,ncol(ordered_compound_percentage_wholemilk)), function(N){
  return(rowSums(ordered_compound_percentage_wholemilk[,1:N]))
})
colnames(sum_top_N_compound_percentage_wholemilk) = c("top_10", "top_20", "top_50", "top_100", "all_compounds")

#Reducedfat sample only
ordered_compound_percentage_reducedfat = ordered_compound_percentage[which(rownames(ordered_compound_percentage) %in% c("DMD100029","DMD100031")),]
sum_top_N_compound_percentage_reducedfat = sapply(c(10,20,50,100,ncol(ordered_compound_percentage_reducedfat)), function(N){
  return(rowSums(ordered_compound_percentage_reducedfat[,1:N]))
})
colnames(sum_top_N_compound_percentage_reducedfat) = c("top_10", "top_20", "top_50", "top_100", "all_compounds")


convert_snp_name = function(snps){
  return(sapply(strsplit(snps,":",fixed=T), function(x){
    partA = sprintf("GK%06d.2",as.numeric(x[1]))
    partB = gsub("/","_",x[2])
    return(paste0(partA,"_",partB))
  }))
}

convert_snp_name2 = function(snps){
  return(sapply(strsplit(snps,"_",fixed=T), function(x){
    partA = gsub("GK","",x[1])
    partA = gsub(".2","",partA,fixed=T)
    partA = as.numeric(partA)
                       
    return(paste0(partA,":",x[2],"_",x[3],"/",x[4]))
  }))
}