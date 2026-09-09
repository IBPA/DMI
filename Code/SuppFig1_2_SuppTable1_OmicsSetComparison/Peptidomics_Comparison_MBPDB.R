#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))

source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

filelist = list.files(paste0(DATA_PATH,"\\DatabaseComparison\\MBPDB"), pattern = "*tsv")
files = lapply(filelist, function(x) {
    read.delim(paste0(DATA_PATH,"\\DatabaseComparison/MBPDB/", x), header = TRUE, stringsAsFactors = FALSE)
  }
)
merged_data = do.call(rbind, lapply(files,function(x){
  x[1:(nrow(x)-1),]
}))

quantified_peptides = sapply(colnames(data_pca_tsne)[which(sapply(colnames(data_pca_tsne), molecular_classification) == "Peptides")],
                             function(x){
                               idx = which(mol_data$`DMD ID` == x)
                               return(mol_data$`Chemical Composition`[idx])
                             })
mbpdb_peptides = merged_data$Peptide

# Compare the two sets of peptides
intersected_peptides = intersect(quantified_peptides, mbpdb_peptides) #134
unique_quantified_peptides = setdiff(quantified_peptides, mbpdb_peptides) #4278
unique_mbpdb_peptides = setdiff(mbpdb_peptides, quantified_peptides) #356

mbpdb_comparison_summary = rbind(cbind(intersected_peptides, rep("In Both Dataset",length(intersected_peptides))),
                                 cbind(unique_quantified_peptides, rep("In Our Dataset Only",length(unique_quantified_peptides))),
                                 cbind(unique_mbpdb_peptides, rep("In MBPDB Only",length(unique_mbpdb_peptides))))

dir.create(paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison"), showWarnings = F)
write.csv(mbpdb_comparison_summary,paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison\\Peptidomics_Comparison_MBPDB.csv"))


#Source protein
#Intersection: a-lactalbumin, a-s1-casein, a-s2-casein, b-casein, k-casein, b-lactoglobulin, 
#hemoglobin subunit alpha, lactotransferrin, glycosylation-depedent cell adhesion molecule 1, serum albumin

#MDPDB only: k-casein genetic variant (4)
#G1
#C
#F1
#B,F2,G1,G2,J
