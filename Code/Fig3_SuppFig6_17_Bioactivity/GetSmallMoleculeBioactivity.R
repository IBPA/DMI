#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

dir.create(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity"), showWarnings = FALSE)

library(jsonlite)
library(PubChemR)
library(xml2)
library(XML)
library(basicdrm)

mol_data_pubchem = lapply(mol_data$`External Database IDs`, function(x){
  if (x == "") return(NA)
  parsed = fromJSON(x)
  if (!is.null(parsed$PubChem)) {
    return(parsed$PubChem)
  } else {
    return(NA)
  }
})
mol_data_pubchem[[2141]] = c("587","5359254","9548602") #Manually fix the typo error for DMD302151 (Creatine Phosphate) => recollecting the data is necessary
names(mol_data_pubchem) = mol_data$`DMD ID`
mol_data_pubchem = mol_data_pubchem[!is.na(mol_data_pubchem)]

cids = unique(unlist(mol_data_pubchem))
cids = cids[!is.na(cids)]
#cid_data_query = NULL
for (i in 1:length(cids)){
  if (i %% 100 == 0) print(i)
  url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",cids[i],"/assaysummary/CSV")
  csv_content = RCurl::getURL(url)
  csv_data = read.csv(text=csv_content, header = T, check.names = F)
  if (nrow(csv_data) == 0) next
  if (ncol(csv_data) == 1 && colnames(csv_data)[1] == "Status: 404") next
  cid_data_query = rbind(cid_data_query, csv_data)
}

cid_data_query = cid_data_query[!duplicated(cid_data_query),]
write.csv(cid_data_query, paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity\\MilkSmallMoleculeBioactivity_All.csv"), row.names = F)  #Pubchem bioassays -- milk small molecule pairs

cid_data_query_active = cid_data_query[cid_data_query$`Activity Outcome` == "Active",]
write.csv(cid_data_query_active, paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity\\MilkSmallMoleculeBioactivity_ActiveOnly.csv"), row.names = F) #Pubchem bioassays -- milk small molecule pairs (Labeled as active only)


aids = unique(cid_data_query_active$AID)
aid_data_query_list = list()
for (i in 1:length(aids)){
  if (i %% 100 == 0) print(i)
  
  cur_cid = cid_data_query_active[cid_data_query_active$AID == aids[i], "CID"]
  if (length(cur_cid < 10)){
    url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/",aids[i],"/CSV?cid=", paste(cur_cid, collapse = ","))
  }else{
    url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/",aids[i],"/CSV")
  }
  
  csv_content = RCurl::getURL(url)
  csv_data = read.csv(text=csv_content, header = T, check.names = F)
  if (nrow(csv_data) == 0) next
  if (ncol(csv_data) == 1 && colnames(csv_data)[1] == "Status: 404") next
  if (ncol(csv_data) == 1) stop('!')
  aid_data_query_list[[i]] = csv_data
}

names(aid_data_query_list) = aids

cid_to_dmd_id = function(cid){
  idx = which(sapply(mol_data_pubchem, function(x) cid %in% x))
  if (length(idx) == 0) {
    print(cid)
    stop('!')
  }
  output_dmd_id = intersect(names(mol_data_pubchem)[idx], colnames(data_pca_tsne))
  if (length(output_dmd_id) == 1) return(output_dmd_id)
  if (length(output_dmd_id) == 0) return(NA) #Remove the compound which are not quantified
  
  
  best_idx = which.max(colMeans(data_pca_tsne[,output_dmd_id,drop=F]))
  return(output_dmd_id[best_idx])
}

cid_to_dmd_id_all = function(cid){
  idx = which(sapply(mol_data_pubchem, function(x) cid %in% x))
  if (length(idx) == 0) {
    print(cid)
    stop('!')
  }
  output_dmd_id = intersect(names(mol_data_pubchem)[idx], colnames(data_pca_tsne))
  if (length(output_dmd_id) == 1) return(output_dmd_id)
  if (length(output_dmd_id) == 0){
    if (length(names(mol_data_pubchem)[idx]) == 1){
      return(names(mol_data_pubchem)[idx])
    }else{
      return(paste0(names(mol_data_pubchem)[idx], collapse = "|"))
    }
  }
  
  return(paste0(output_dmd_id, collapse = "|"))
}


#Data with Fitted Hill Curve
aid_data_fitted_curve_list = list()
name_aid_fitted_curve = c()
target_column = c("Fit_LogAC50","Fit_HillSlope", "Fit_R2", "Fit_InfiniteActivity","Fit_ZeroActivity","Fit_CurveClass")
for (i in 1 : length(aid_data_query_list)){
  if (i %% 100 == 0) print(i)
  if (any(grepl("Fit_",colnames(aid_data_query_list[[i]])))){
    basic_info = aid_data_query_list[[i]][,grepl("PUBCHEM_",colnames(aid_data_query_list[[i]]))]
    raw_fit_info = aid_data_query_list[[i]][,grepl("Fit_",colnames(aid_data_query_list[[i]]))]
    efficacy_info = aid_data_query_list[[i]][,grepl("Efficacy",colnames(aid_data_query_list[[i]])),drop=F]
    if (ncol(efficacy_info) == 0) next
    
    if (length(setdiff(target_column, colnames(raw_fit_info))) == 0){
      final_fit_info = raw_fit_info[,target_column]
      efficacy_info = aid_data_query_list[[i]][,"Efficacy",drop=F]
      
      new_data = cbind(basic_info, final_fit_info, efficacy_info)
      new_data = new_data[which(new_data[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Active"),,drop=F]
      new_data = new_data[which(new_data[,"PUBCHEM_CID"] %in% cids),,drop=F]
      new_data = cbind(names(aid_data_query_list)[i], NA, sapply(new_data[,"PUBCHEM_CID"], cid_to_dmd_id), new_data)
      colnames(new_data)[1:3] = c("AID", "Rep", "DMD_ID")
      
      new_data = new_data[!is.na(new_data[,"DMD_ID"]),,drop=F]
      new_data = new_data[!apply(new_data[,target_column],1,function(x){any(x=="")}),,drop=F]
      if (nrow(new_data) > 0){
        aid_data_fitted_curve_list[[length(aid_data_fitted_curve_list)+1]] = new_data
        name_aid_fitted_curve = c(name_aid_fitted_curve, names(aid_data_query_list)[i])
      }
    }else{
      rep_idx = 1
      while(TRUE){
        cur_target_column = paste0(target_column,"-Replicate_",rep_idx)
        cur_efficacy_column = paste0("Efficacy-Replicate_",rep_idx)
        
        if (length(setdiff(cur_target_column, colnames(raw_fit_info))) == 0){
          final_fit_info = raw_fit_info[,cur_target_column]
          colnames(final_fit_info) = target_column
          
          efficacy_info = aid_data_query_list[[i]][,cur_efficacy_column,drop=F]
          colnames(efficacy_info) = "Efficacy"
          
          new_data = cbind(basic_info, final_fit_info, efficacy_info)
          new_data = new_data[which(new_data[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Active"),,drop=F]
          new_data = new_data[which(new_data[,"PUBCHEM_CID"] %in% cids),,drop=F]
          new_data = cbind(names(aid_data_query_list)[i], rep_idx, sapply(new_data[,"PUBCHEM_CID"], cid_to_dmd_id), new_data)
          colnames(new_data)[1:3] = c("AID", "Rep", "DMD_ID")
          
          new_data = new_data[!is.na(new_data[,"DMD_ID"]),,drop=F]
          new_data = new_data[!apply(new_data[,target_column],1,function(x){any(x=="")}),,drop=F]
          
          if (nrow(new_data) > 0){
            aid_data_fitted_curve_list[[length(aid_data_fitted_curve_list)+1]] = new_data
            name_aid_fitted_curve = c(name_aid_fitted_curve, paste0(names(aid_data_query_list)[i],"-Rep_",rep_idx))
          }
        }else{
          break
        }
        rep_idx = rep_idx + 1
      }
    }
  }
}
names(aid_data_fitted_curve_list) = name_aid_fitted_curve


merged_aid_data_fitted_curve = do.call(rbind, aid_data_fitted_curve_list)
write.csv(merged_aid_data_fitted_curve, paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity\\MilkSmallMoleculeBioactivity_FittedCurve.csv"), row.names = F) #Pubchem bioassays -- milk small molecule pairs (Labeled as active and contain fitted Hill Curve only)


#Extend the bioassays with Standard Type (AC50, EC50, IC50) -- the efficacy can be calculated based on the representative Hill Curve (Supp. Fig 17)
aid_data_standard_type_list = list()
name_aid_standard_type = c()
target_column = c("PubChem Standard Value", "Standard Type")
for (i in 1 : length(aid_data_query_list)){
  if (i %% 100 == 0) print(i)
  basic_info = aid_data_query_list[[i]][,grepl("PUBCHEM_",colnames(aid_data_query_list[[i]]))]
  
  if (length(setdiff(target_column, colnames(aid_data_query_list[[i]]))) == 0){
    standard_info = aid_data_query_list[[i]][,target_column]
    
    new_data = cbind(basic_info, standard_info)
    new_data = new_data[which(new_data[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Active"),,drop=F]
    new_data = new_data[which(new_data[,"PUBCHEM_CID"] %in% cids),,drop=F]
    new_data = cbind(names(aid_data_query_list)[i], NA, sapply(new_data[,"PUBCHEM_CID"], cid_to_dmd_id), new_data)
    colnames(new_data)[1:3] = c("AID", "Rep", "DMD_ID")
    
    new_data = new_data[!is.na(new_data[,"DMD_ID"]),,drop=F]
    new_data = new_data[!apply(new_data[,target_column],1,function(x){any(x=="")}),,drop=F]
    if (nrow(new_data) > 0){
      aid_data_standard_type_list[[length(aid_data_standard_type_list)+1]] = new_data
      name_aid_standard_type = c(name_aid_standard_type, names(aid_data_query_list)[i])
    }
  }else{
    rep_idx = 1
    while(TRUE){
      cur_target_column = paste0(target_column,"-Replicate_",rep_idx)
      
      if (length(setdiff(cur_target_column, colnames(aid_data_query_list[[i]]))) == 0){
        standard_info = aid_data_query_list[[i]][,cur_target_column]
        colnames(final_fit_info) = target_column
        
        new_data = cbind(basic_info, standard_info)
        new_data = new_data[which(new_data[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Active"),,drop=F]
        new_data = new_data[which(new_data[,"PUBCHEM_CID"] %in% cids),,drop=F]
        new_data = cbind(names(aid_data_query_list)[i], rep_idx, sapply(new_data[,"PUBCHEM_CID"], cid_to_dmd_id), new_data)
        colnames(new_data)[1:3] = c("AID", "Rep", "DMD_ID")
        
        new_data = new_data[!is.na(new_data[,"DMD_ID"]),,drop=F]
        new_data = new_data[!apply(new_data[,target_column],1,function(x){any(x=="")}),,drop=F]
        
        if (nrow(new_data) > 0){
          aid_data_standard_type_list[[length(aid_data_standard_type_list)+1]] = new_data
          name_aid_standard_type = c(name_aid_standard_type, paste0(names(aid_data_query_list)[i],"-Rep_",rep_idx))
        }
      }else{
        break
      }
      rep_idx = rep_idx + 1
    }
  }
}
names(aid_data_standard_type_list) = name_aid_standard_type


merged_aid_data_standard_type = do.call(rbind, aid_data_standard_type_list)
merged_aid_data_standard_type_XC50 = merged_aid_data_standard_type[which(merged_aid_data_standard_type[,"Standard Type"] %in% c("AC50","EC50","IC50")),,drop=F]
write.csv(merged_aid_data_standard_type_XC50, paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity\\MilkSmallMoleculeBioactivity_StandardType.csv"), row.names = F) #The bioassays (contains milk small molecules) with standard type of AC50, EC50, IC50 (Labeled as active and contain standard type only)


#Get the bioassay description for selecting the health/disease relevant bioassays
#These bioassays should contains the fitted Hill curve or should be measured under the standard type (AC50, EC50, IC50)
aids_xc50_curve = union(unique(merged_aid_data_standard_type_XC50$AID),unique(merged_aid_data_fitted_curve$AID))
aid_desc_posctrl = c()
aid_desc_ctrl = c()
aid_desc_normalize = c()
name_aid_ctrl = c()

pattern_posctrl = "[Pp]ositive [Cc]ontrol"
pattern_ctrl = "[Cc]ontrol"
pattern_normalize = "[Nn]ormalize"

for (i in 1 : length(aids_xc50_curve)){
  if (i %% 100 == 0) print(i)
  url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/",aids_xc50_curve[i],"/description/XML")
  xml_content = RCurl::getURL(url)
  xml_data = read_xml(xml_content)
  
  xml_text = xml_text(xml_data)
  
  if (grepl(pattern_normalize, xml_text) || grepl(pattern_ctrl, xml_text)){
    if (grepl(pattern_normalize, xml_text)){
      aid_desc_normalize = c(aid_desc_normalize, str_extract(xml_text, paste0(".{0,300}(?=",pattern_normalize,").{0,300}")))
    }else{
      aid_desc_normalize = c(aid_desc_normalize, NA)
    }
    if (grepl(pattern_ctrl, xml_text(xml_data))){
      if (grepl(pattern_posctrl, xml_text(xml_data))){
        aid_desc_posctrl = c(aid_desc_posctrl, str_extract(xml_text(xml_data), paste0(".{0,300}(?=",pattern_posctrl,").{0,300}")))
        aid_desc_ctrl = c(aid_desc_ctrl, str_extract(xml_text(xml_data), paste0(".{0,300}(?=",pattern_posctrl,").{0,300}")))
      }else{
        aid_desc_posctrl = c(aid_desc_posctrl, NA)
        aid_desc_ctrl = c(aid_desc_ctrl, str_extract(xml_text(xml_data), paste0(".{0,300}(?=",pattern_ctrl,").{0,300}")))
      }
    }else{
      aid_desc_posctrl = c(aid_desc_posctrl, NA)
      aid_desc_ctrl = c(aid_desc_ctrl, NA)
    }
    name_aid_ctrl = c(name_aid_ctrl, aids_xc50_curve[i])
  }
  
}
write.csv(data.frame(AID = name_aid_ctrl, Description_Normalize = aid_desc_normalize, Description_PosCtrl = aid_desc_posctrl, Description_Ctrl = aid_desc_ctrl), 
          paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_Bioactivity\\BioassayDescription.csv", row.names = F))


#The Positive Control, Negative Control, Definition of 0%/100% efficacy, and purpose of the bioassays were then manually annotated (See the file: BioassayDescription_Manual.csv as reference)
#The health/disease relevant bioassays will be selected manually (See the file: BioassayDescription_Manual_Selected.csv as reference)



