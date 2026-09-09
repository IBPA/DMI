#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))

source(paste0(CODE_PATH,"\\Analysis_Utils.r"))
library(Biostrings)

#Load BoMiProt data
fasta_files = list.files(paste0(DATA_PATH,"\\DatabaseComparison\\BoMiProt\\bomiprot_fasta"), pattern = "\\.fasta$", full.names = T)


bomiprot_seq = sapply(fasta_files, function(fasta_file){
    protein_seqs = readAAStringSet(fasta_file)
    cur_protein_seq = as.character(protein_seqs[[1]])
    return (cur_protein_seq)
})

mol_data_protein = mol_data[which(grepl(mol_data$`Omic Lab`, pattern = "Proteomics")),]
#The protein not in BoMiProt (Supp Fig 1)
mol_data_protein_new = mol_data_protein[which(!(mol_data_protein$`Chemical Composition` %in% bomiprot_seq)),]

intersect_protein = mol_data_protein$`Chemical Composition`[which((mol_data_protein$`Chemical Composition` %in% bomiprot_seq))]
not_in_BoMiProt = mol_data_protein$`Chemical Composition`[which((!mol_data_protein$`Chemical Composition` %in% bomiprot_seq))]
BoMiProtOnly = setdiff(bomiprot_seq, mol_data_protein$`Chemical Composition`)




#The proteins in the top 100 most abundant quantified compound set and not in BoMiProt (Supp Table 1)
dmd_id_top100_new = intersect(colnames(ordered_compound_percentage)[1:100], mol_data_protein_new$`DMD ID`)
mol_data_protein_top100_new = lapply(dmd_id_top100_new, function(x){
    idx = which(mol_data_protein_new$`DMD ID` == x)
    return(mol_data_protein_new[idx,])
})
mol_data_protein_top100_new = do.call(rbind, mol_data_protein_top100_new)

rank = match(dmd_id_top100_new, colnames(ordered_compound_percentage))
avg_conc = apply(data_pca_tsne[,dmd_id_top100_new], 2, mean)/240
avg_content = apply(ordered_compound_percentage[,dmd_id_top100_new], 2, mean)*100

min_conc = apply(data_pca_tsne[,dmd_id_top100_new], 2, min)/240
max_conc = apply(data_pca_tsne[,dmd_id_top100_new], 2, max)/240
min_content = apply(ordered_compound_percentage[,dmd_id_top100_new], 2, min)*100
max_content = apply(ordered_compound_percentage[,dmd_id_top100_new], 2, max)*100

#Here is the summary of the proteins in the top 100 most abundant quantified compound set and not in BoMiProt
summary_protein_top100_new = cbind(mol_data_protein_top100_new, rank, avg_conc, avg_content, min_conc, max_conc, min_content, max_content)

dir.create(paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison"), showWarnings = F)
BoMiProt_comparison_summary = rbind(cbind(intersect_protein, rep("In Both Dataset",length(intersect_protein))),
                                    cbind(not_in_BoMiProt, rep("In Our Dataset Only",length(not_in_BoMiProt))),
                                    cbind(BoMiProtOnly, rep("In BoMiProt Only",length(BoMiProtOnly))))


write.csv(BoMiProt_comparison_summary,paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison\\Proteomics_Comparison_BoMiProt.csv"))
write.csv(summary_protein_top100_new,paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison\\Protein_Top100_NotInBoMiProt.csv"))



BoMiProt_comparison_summary
summary_protein_top100_new