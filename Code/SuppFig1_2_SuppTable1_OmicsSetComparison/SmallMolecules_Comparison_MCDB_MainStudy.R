#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))

source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

library(gsubfn)
library(xml2)
xml_data = read_xml(paste0(DATA_PATH,"\\DatabaseComparison/MCDB/milk_metabolites/milk_metabolites.xml"))
metabolite_nodes = xml_children(xml_data)

accessions = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  accession_idx = which(sapply(children, function(child) xml_name(child) == "accession"))
  if (length(accession_idx) != 1) stop()
  accession = xml_text(children[accession_idx])
  return(accession)
} )
inchikeys = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  inchikey_idx = which(sapply(children, function(child) xml_name(child) == "inchikey"))
  if (length(inchikey_idx) != 1) stop()
  inchikey = xml_text(children[inchikey_idx])
  return(inchikey)
} )
smiles = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  smile_idx = which(sapply(children, function(child) xml_name(child) == "smiles"))
  if (length(smile_idx) != 1) stop()
  smile = xml_text(children[smile_idx])
  return(smile)
} )
pubchem_cids = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  pubchem_cid_idx = which(sapply(children, function(child) xml_name(child) == "pubchem_compound_id"))
  if (length(pubchem_cid_idx) == 0) return(NA)
  if (length(pubchem_cid_idx) != 1) stop()
  pubchem_cid = xml_text(children[pubchem_cid_idx])
  return(pubchem_cid)
} )
cas = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  cas_idx = which(sapply(children, function(child) xml_name(child) == "cas_registry_number"))
  if (length(cas_idx) != 1) stop()
  cas = xml_text(children[cas_idx])
  return(cas)
} )

names = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  name_idx = which(sapply(children, function(child) xml_name(child) == "name"))
  if (length(name_idx) != 1) stop()
  name = xml_text(children[name_idx])
  return(name)
} )

synonyms = lapply(metabolite_nodes, function(node){
  children = xml_children(node)
  synonym_idx = which(sapply(children, function(child) xml_name(child) == "synonyms"))
  if (length(synonym_idx) != 1) stop()
  synonym = xml_text(xml_children(children[synonym_idx]))
  return(synonym)
} )

concentrations_ind = sapply(metabolite_nodes, function(node){
  children = xml_children(node)
  conc_idx = which(sapply(children, function(child) xml_name(child) == "concentrations"))
  if (length(conc_idx) > 1) stop()
  if (length(conc_idx) == 0) return("No record in the main study")
  
  found_in_main_study = F
  
  concentration_fields = xml_children(children[conc_idx])
  for (cur_concentration_field in concentration_fields){
    cur_found_in_main_study = F
    
    conc_reference_idx = which(sapply(xml_children(cur_concentration_field ), function(child) xml_name(child) == "references"))
    if (length(conc_reference_idx) > 0){
      available_references = xml_children((xml_children(cur_concentration_field)[conc_reference_idx]))
      for (cur_available_references in available_references){
        title = xml_text(xml_children(cur_available_references)[1])
        if (length(xml_children(cur_available_references)) == 2){
          pubmed_id = xml_text(xml_children(cur_available_references)[2])
        }else{
          pubmed_id = ""
        }
        
        
        if (title == "A. Foroutan et al. The Chemical Composition of Cow’s Milk (in preparation)" || pubmed_id == "30994344"){
          found_in_main_study = T
          cur_found_in_main_study = T
          break
        }
      }
    }
    
    
    conc_value_idx = which(sapply(xml_children(cur_concentration_field ), function(child) xml_name(child) == "concentration_value"))
    conc_units_idx = which(sapply(xml_children(cur_concentration_field ), function(child) xml_name(child) == "concentration_units"))
    
    
    if (length(conc_value_idx) == 0 || length(conc_units_idx) == 0) next
    if (length(conc_value_idx) > 1 || length(conc_units_idx) > 1) stop()
    conc_value = xml_text(xml_children(cur_concentration_field)[conc_value_idx])
    conc_units = xml_text(xml_children(cur_concentration_field)[conc_units_idx])
    if (conc_value == "" || conc_units == "") next
    
    if (cur_found_in_main_study){
      print(paste0("Found concentration in the main study: ", conc_value, " ", conc_units))
      return("Found Concentration in Main Study")
    }
  }
  
  if (found_in_main_study){
    return("Identified in Main Study")
  }else{
    return("No record in the main study")
  }
  
} )

lower_mcdb_synonyms = lapply(synonyms, function(syn_list){
  if (length(syn_list) == 0) return("")
  syn_list = gsubfn(" ","-", syn_list, fixed=T)
  return(gsubfn("[A-Z]", tolower, syn_list))
})
lower_names = gsubfn("[A-Z]", tolower, names)
lower_names = gsubfn(" ","-", lower_names, fixed=T)
lower_synonyms = lapply(1 : nrow(mol_data), function(i){
  name = mol_data$`Molecule Name`[i]
  display_name = mol_data$`Display Name`[i]
  synonyms = gsubfn("{","",mol_data$Synonyms[i], fixed=T)
  synonyms = gsubfn("}","",synonyms, fixed=T)
  if (synonyms != ""){
    synonyms = as.character(read.csv(text=synonyms, header = F, stringsAsFactors = F))
    synonyms = gsub(" ","-",synonyms, fixed=T)
  }
  final_name_list = c(name, display_name, synonyms)
  return(gsubfn("[A-Z]", tolower, final_name_list))
})
names(lower_synonyms) = mol_data$`DMD ID`


direct_match_attempt = sapply(1:length(lower_names), function(i){
  name = lower_names[i]
  mcdb_synonyms = lower_mcdb_synonyms[[i]]
  match_ind = sapply(lower_synonyms, function(syn_list){
    return((name %in% syn_list && name != "") || length(setdiff(intersect(mcdb_synonyms, syn_list),"")) > 0)
  })
  
  idx = which(match_ind == T)
  if (length(idx) == 0) return(NA)
  if (length(idx) > 1) return(paste(mol_data$`DMD ID`[idx], collapse=";"))
  
  return(mol_data$`DMD ID`[idx])
})

mcdb_summary = cbind(accessions, names, inchikeys, smiles, pubchem_cids, concentrations_ind, direct_match_attempt)
colnames(mcdb_summary)[7] = "DMD_ID_Direct_Name_Match"

mol_data_small_molecules_detected = mol_data[which(sapply(mol_data$`DMD ID`, molecular_classification) %in% c("Lipids","Metabolites","Oligosaccharides","Metabolites (Untargeted)")),]
mol_data_small_molecules_detected_pubchem = lapply(mol_data_small_molecules_detected$`External Database IDs`, function(x){
  if (x == "") return(NA)
  parsed = fromJSON(x)
  if (!is.null(parsed$PubChem)) {
    return(parsed$PubChem)
  } else {
    return(NA)
  }
})
names(mol_data_small_molecules_detected_pubchem) = mol_data_small_molecules_detected$`DMD ID`
mol_data_small_molecules_detected_pubchem[['DMD302151']] = c("587","5359254","9548602")


find_mapping_dmd_id = t(sapply(1 : nrow(mcdb_summary), function(i){
  cur_pubchem_id = mcdb_summary[i,"pubchem_cids"]
  invalid_pubchem_id = F
  if (is.na(cur_pubchem_id) || (cur_pubchem_id == "")){
    invalid_pubchem_id = T
  }else{
    invalid_pubchem_id = F
  }
  idx = which(sapply(mol_data_small_molecules_detected_pubchem, function(x) cur_pubchem_id %in% x))
  if (length(idx) == 0 || invalid_pubchem_id) {
    cur_name_mapping_ids = strsplit(mcdb_summary[i,"DMD_ID_Direct_Name_Match"], ";")[[1]]
    if (is.na(cur_name_mapping_ids[1])){
      return(c(NA, "No Match"))
    }
    if (length(cur_name_mapping_ids) == 1){
      return(c(cur_name_mapping_ids, "Unique Match (by names only)"))
    }
    return(c(as.character(mcdb_summary[i,"DMD_ID_Direct_Name_Match"]), "Multiple Match (by names only)"))
  }else{
    cur_name_mapping_ids = strsplit(mcdb_summary[i,"DMD_ID_Direct_Name_Match"], ";")[[1]]
    cur_ids = names(mol_data_small_molecules_detected_pubchem)[idx]
    
    if (length(intersect(cur_name_mapping_ids, cur_ids)) > 0){
      if (length(intersect(cur_name_mapping_ids, cur_ids)) == 1){
        return(c(intersect(cur_name_mapping_ids, cur_ids), "Unique Match (by names and pubchem)"))
      }
      return(c(paste(intersect(cur_name_mapping_ids, cur_ids), collapse=';'), "Multiple Match (by names and pubchem)"))
    }else{
      if (length(cur_ids) == 1){
        return(c(cur_ids, "Unique Match (by pubchem only)"))
      }
      return(c(paste(cur_ids, collapse=";"), "Multiple Match (by pubchem only)"))
    }
    
  }
  #Should not happen
}))


mcdb_comparison_result = cbind(mcdb_summary, find_mapping_dmd_id)
mcdb_comparison_result

#The following steps is for counting the 
all_mapped_compounds_in_dmd = setdiff(union(setdiff(unique(unlist(sapply(mcdb_data$Final,function(x){strsplit(x,";",fixed=T)}))),NA), unique(unlist(sapply(mcdb_data$DMD_ID_Direct_Name_Match,function(x){strsplit(x,";",fixed=T)})))), NA)
dmd_only = setdiff(mol_data$`DMD ID`[which(grepl("Targeted Metabolomics",mol_data$`Omic Lab`,fixed=T))], all_mapped_compounds_in_dmd)

dir.create(paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison"), showWarnings = F)
write.csv(mcdb_comparison_result,paste0(PROCESSED_DATA_PATH,"\\SuppFig1_2_SuppTable1_OmicsSetComparison\\SmallMolecules_Comparison_MCDB_MainStudy.csv"))

#Note: The comparison results between our study and MCDB are double-checked manually.
#      The MCDB also performed their own study and it is the largest small molecule profiling dataset.
#      The comparison results of this main study is also included in the merged, manually checked table.
#(The manually checked result is provided in the same directory (SmallMolecules_Comparison_MCDB_merged_manual.csv)).


