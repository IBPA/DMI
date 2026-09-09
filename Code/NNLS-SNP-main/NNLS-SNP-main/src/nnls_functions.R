run_nnls_merge = function(alt_ratio, ref_breed_matrix){
  library(nnls)
  alt_ratio = alt_ratio[!is.na(alt_ratio)]
  alt_ratio_entire = rep(mean(alt_ratio),nrow(ref_breed_matrix))
  names(alt_ratio_entire) = rownames(ref_breed_matrix)
  intersect_snp = intersect(rownames(ref_breed_matrix),names(alt_ratio))
  alt_ratio_entire[intersect_snp] = alt_ratio[intersect_snp]
  
  ref_breed_matrix_intersect = ref_breed_matrix[intersect_snp,]
  alt_ratio_intersect = alt_ratio[intersect_snp]
  
  nnls_result_intersect = nnls(rbind(ref_breed_matrix_intersect,
                                     rep(length(alt_ratio_intersect),ncol(ref_breed_matrix))),
                               c(alt_ratio_intersect, 
                                 length(alt_ratio_intersect)))
  
  nnls_result_list = list()
  nnls_result_list[['Intersect']] = nnls_result_intersect
  nnls_result_list[['AltRatioFilled']] = alt_ratio_entire
  nnls_result_list[['IntersectSNP']] = intersect_snp
  
  cor_intersect = cor(ref_breed_matrix_intersect, alt_ratio_intersect)
  nnls_result_list[['cor_intersect']] = cor_intersect
  
  return(nnls_result_list)
}



run_nnls_merge_batch = function(ad_matrix, ref_matrix){
  results = lapply(1:ncol(ad_matrix), function(i){
    result_this_sample = run_nnls_merge(ad_matrix[,i], ref_matrix)
    return(result_this_sample$Intersect$x)
  })
  
  results = do.call(rbind,results)
  
  rownames(results) = colnames(ad_matrix)
  colnames(results) = colnames(ref_matrix)
  
  return(results)
}

library(foreach)
library(doParallel)
registerDoParallel(4)

subsample_run_nnls_merge_batch = function(ad_matrix, ref_matrix, n_rep = 100, n_snp = c(10, 50, 1e2, 1e3, 1e4, 1e5, 1e6, Inf)){
  result_list = list()
  n_snp[is.infinite(n_snp)] = nrow(ad_matrix)
  result_list[['N_SNP']] = n_snp
  result_list[['n_rep']] = n_rep
  result_list[['nnls_result']] = list()
  for (i in 1 : length(n_snp)){
    print(paste0("Subsampling: Current #SNP = ", n_snp[i]))
    result_list[['nnls_result']][[i]] = list()
    if (n_snp[i] == nrow(ad_matrix)){
      result_list[['nnls_result']][[i]][[1]] = run_nnls_merge_batch(ad_matrix, ref_matrix)
    }else{
      result_list[['nnls_result']][[i]] = foreach (j = 1 : n_rep) %dopar% {
        sample_idx = sample(1:nrow(ad_matrix),n_snp[i])
        run_nnls_merge_batch(ad_matrix[sample_idx,], ref_matrix[sample_idx,])
      }
    }
    
  }
  return(result_list)
}

extend_sample_run_nnls_merge_batch = function(ad_matrix, ref_matrix, ad_matrix_extend, ref_matrix_extend, 
                                              n_rep = 100, n_snp = c(10, 50, 1e2, 1e3, 1e4, 1e5, 1e6)){
  result_list = list()
  result_list[['N_SNP']] = n_snp
  result_list[['n_rep']] = n_rep
  result_list[['nnls_result']] = list()
  
  addl_snps = setdiff(rownames(ad_matrix_extend), rownames(ad_matrix))
  ad_matrix_addl_snps = ad_matrix_extend[addl_snps,]
  ref_matrix_addl_snps = ref_matrix_extend[addl_snps,]
  
  for (i in 1 : length(n_snp)){
    print(paste0("Subsampling: Current #SNP = ", n_snp[i]))
    if (n_snp[i] <= nrow(ad_matrix)){
      result_list[['nnls_result']][[i]] = foreach (j = 1 : n_rep) %dopar% {
        sample_idx = sample(1:nrow(ad_matrix),n_snp[i])
        run_nnls_merge_batch(ad_matrix[sample_idx,], ref_matrix[sample_idx,])
      }
    }else{
      result_list[['nnls_result']][[i]] = foreach (j = 1 : n_rep) %dopar% {
        n_addl_snps = n_snp[i] - nrow(ad_matrix)
        sample_idx = sample(1:nrow(ad_matrix_addl_snps),n_addl_snps)
        
        ad_matrix_in = rbind(ad_matrix, ad_matrix_addl_snps[sample_idx,])
        ref_matrix_in = rbind(ref_matrix, ref_matrix_addl_snps[sample_idx,])
        
        run_nnls_merge_batch(ad_matrix_in, ref_matrix_in)
      }
    }
  }
  return(result_list)
}



extract_ad_matrix_info_batch = function(filepath){
  batch_data = read.vcfR(filepath)
  batch_data_ad = extract.gt(batch_data, element="AD")
  chrom_new_label = mapping_dict[batch_data@fix[,"CHROM"]]
  na_ind_new_label = is.na(chrom_new_label)
  chrom_new_label[na_ind_new_label] = batch_data@fix[na_ind_new_label,"CHROM"]
  chrom_new_label = as.character(chrom_new_label)
  
  entry_tag = paste0(chrom_new_label,":",batch_data@fix[,"POS"],"_",
                     batch_data@fix[,"REF"],"/",batch_data@fix[,"ALT"])
  
  tmp_split = lapply(1:ncol(batch_data_ad),function(i){strsplit(batch_data_ad[,i],',',fixed=T)})
  genotype_count = sapply(1:ncol(batch_data_ad),function(i){sapply(tmp_split[[i]],length)})
  ind_two_genotype_only = apply(genotype_count,1,max) == 2
  
  batch_data_ad_matrix_ref = sapply(tmp_split,function(x){sapply(x,function(y){y[1]})[ind_two_genotype_only]})
  batch_data_ad_matrix_alt = sapply(tmp_split,function(x){sapply(x,function(y){y[2]})[ind_two_genotype_only]})
  
  batch_data_ad_matrix_ref = apply(batch_data_ad_matrix_ref, 2, as.numeric)
  batch_data_ad_matrix_alt = apply(batch_data_ad_matrix_alt, 2, as.numeric)
  
  rownames(batch_data_ad_matrix_ref) = entry_tag[ind_two_genotype_only]
  colnames(batch_data_ad_matrix_ref) = colnames(batch_data_ad)
  
  rownames(batch_data_ad_matrix_alt) = entry_tag[ind_two_genotype_only]
  colnames(batch_data_ad_matrix_alt) = colnames(batch_data_ad)
  
  result = list()
  result[['raw_ad']] = batch_data_ad
  result[['ad_matrix_ref']] = batch_data_ad_matrix_ref
  result[['ad_matrix_alt']] = batch_data_ad_matrix_alt
  result[['ad_matrix']] = batch_data_ad_matrix_alt/(batch_data_ad_matrix_ref + batch_data_ad_matrix_alt)
  
  return(result)
}