#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggrastr)

if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\SuppFig18To23_MOFA\\ad_matrix_snp_medium.rdata"))){
  stop("Missing required allele frequency matrix of the selected SNPs (ad_matrix_snp_medium.rdata).")
}
load(paste0(PROCESSED_DATA_PATH,"\\SuppFig18To23_MOFA\\ad_matrix_snp_medium.rdata"))


snps_significant_medium = rownames(ad_matrix_snp_medium)
pca_compound = prcomp(ordered_compound_percentage, scale. = F, center=T)

#Find the SNPs that are most correlated with the compound PCs
cor_medium_snp_pcs = apply(ad_matrix_snp_medium[,rownames(pca_compound$x)], 1, function(x){
  return(sapply(1:nrow(ordered_compound_percentage),function(y){
    cor(x, pca_compound$x[,y])
  }))
})
cor_medium_snp_pcs[is.na(cor_medium_snp_pcs)] = 0

spearman_medium_snp_pcs = apply(ad_matrix_snp_medium[,rownames(pca_compound$x)], 1, function(x){
  return(sapply(1:nrow(ordered_compound_percentage),function(y){
    cor(x, pca_compound$x[,y], method = "spearman")
  }))
})
spearman_medium_snp_pcs[is.na(spearman_medium_snp_pcs)] = 0

ind_045 = (abs(cor_medium_snp_pcs) > 0.45) & (abs(spearman_medium_snp_pcs) > 0.45)
n_pcs = 7
selected_snps_idx = apply(ind_045[1:n_pcs,,drop=F], 1, function(x){which(x==T)}) #61 SNPs

#PC prediction (61 SNPs -> 7 PCs)
cur_linear_models = lapply(1:n_pcs, function(cur_pc){
  lm(pca_compound$x[,cur_pc] ~ ., as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[cur_pc]]),rownames(pca_compound$x),drop=F])))
})
predicted_pcs = lapply(1:n_pcs, function(cur_pc){
  predict.lm(cur_linear_models[[cur_pc]], as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[cur_pc]]),rownames(pca_compound$x),drop=F])))
})

#Plot PCs (predicted vs actual)
for (i in 1 : n_pcs){
  svg(paste0(FIGURE_PATH,"\\SuppFig25_PC",i,"_SNP_Predict.svg"), width = 3, height = 3)
  plot_df_predicted = data.frame(Actual = pca_compound$x[,i], Predicted = predicted_pcs[[i]], 
                                 Region = factor(regions[product_data$`Region Code`[match(rownames(pca_compound$x), product_data$`DMD ID`)]],levels=regions),
                                 Brand = c("Brand A","Brand B")[product_data$`Brand Code`[match(rownames(pca_compound$x), product_data$`DMD ID`)]],
                                 Time = factor(c("T1","T2","T3")[product_data$`Time Point Code`[match(rownames(pca_compound$x), product_data$`DMD ID`)]], levels=c("T1","T2","T3")))
  p1 = ggplot(plot_df_predicted, aes(x=Actual, y=Predicted, alpha=Time, shape = Brand, fill = Region)) +
    geom_point(size=3, stroke=1.0) +
    theme_classic() +
    labs(title=paste0("PC",i," (",round(summary(pca_compound)$importance[2,i]*100,3),"% variance explained)"), x="Actual PC", y="Predicted PC") +
    theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
          axis.title = element_text(size=10), axis.text = element_text(size=8)) + 
    scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) + 
    scale_shape_manual(values = c(21, 24)) + 
    scale_alpha_manual(values = c(0.3, 0.65, 1.0)) + 
    guides(fill=guide_legend(override.aes=list(shape=21)), 
           alpha=guide_legend(override.aes=list(shape=21, fill="black")), 
           shape=guide_legend(override.aes=list(alpha=1))) + 
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") + 
    geom_text(label=paste0("PCC = ", round(cor(plot_df_predicted$Actual, plot_df_predicted$Predicted), 3), "\nSpearman = ", round(cor(plot_df_predicted$Actual, plot_df_predicted$Predicted, method="spearman"), 3)),
              x = max(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)), 
              y = min(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)), hjust = 1, vjust = 0,inherit.aes = F, check_overlap = T, size=3) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(min(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)), max(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(min(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)), max(c(plot_df_predicted$Actual,plot_df_predicted$Predicted)))) 
  
  print(p1)
  dev.off()
}

#Plot Linear Model Loadings (including intercept)
for (i in 1 : n_pcs){
  svg(paste0(FIGURE_PATH,"\\SuppFig26_PC",i,"_SNP_Predict_Loadings.svg"), width = 3, height = 3)
  loadings = summary(cur_linear_models[[i]])$coefficients[,1]
  names(loadings) = sapply(strsplit(names(loadings),"_",fixed=T), function(x){
    if (length(x) == 1) return("Intercept")
    part_A = gsub("GK","",x[1],fixed=T)
    part_A = as.numeric(gsub(".2","",part_A,fixed=T))
    return(paste0(part_A,":",x[2],"_",x[3],"/",x[4]))
  })
  loadings_df = data.frame(SNP = names(loadings), Loading = loadings)
  loadings_df$SNP = factor(loadings_df$SNP, levels=loadings_df$SNP[order(abs(loadings_df$Loading), decreasing = F)])
  loadings_df$hjust = ifelse(loadings_df$Loading[1] > 0, -1, 1)
  #loadings_df$hjust[1] = 0.5
  
  loadings_df$type = c("Intercept", rep("Positive SNP", nrow(loadings_df)-1))
  loadings_df$type[which(loadings_df$Loading < 0 & loadings_df$type == "Positive SNP")] = "Negative SNP"
  loadings_df$type = factor(loadings_df$type, levels=c("Positive SNP", "Negative SNP", "Intercept"))
  
  loadings_df$Label = paste0(loadings_df$SNP, " : ", round(loadings_df$Loading, 4))
  
  loadings_df$LabelPos = loadings_df$Loading[1]
  loadings_df$LabelPos[1] = loadings_df$Loading[1]/2
  
  loadings_df$Weight = ifelse(abs(loadings_df$Loading) > 1e-4, round(loadings_df$Loading, 4), round(loadings_df$Loading, 6))
  loadings_df$WeightPos = ifelse(loadings_df$Loading > 0, loadings_df$Loading + max(abs(loadings_df$Loading))*0.5, loadings_df$Loading - max(abs(loadings_df$Loading))*0.5)
  
  range_loadings = max(loadings_df$Loading) - min(loadings_df$Loading)
  p1 = ggplot(loadings_df, aes(y=SNP, x=Loading)) +
    geom_bar(stat="identity", aes(fill=type)) +
    scale_fill_manual(values = c("Positive SNP" = "steelblue", "Negative SNP" = "#E41A1C", "Intercept" = "grey")) +
    theme_classic() +
    labs(title="", y="", x="") +
    theme(plot.title = element_text(hjust = 0.5, size=10), 
          #axis.text.x = element_blank(), 
          #axis.ticks.x = element_blank(), 
          axis.title = element_text(size=10), 
          axis.text.y = element_text(size=8),
          legend.position = "none") + 
    geom_hline(yintercept = 0, linetype="dashed", color="red") + 
    #geom_text(aes(label=Label, x=LabelPos, y=SNP, hjust=hjust), vjust=0.5, size=2.5) + 
    geom_text(aes(label=Weight, x=WeightPos, y=SNP), vjust=0.5, size=2.5) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(min(loadings_df$Loading)-range_loadings*0.6, max(loadings_df$Loading)+range_loadings*0.6)) +
    geom_vline(xintercept = 0, color="black")
  print(p1)
  dev.off()
}

#Performance Evaluation
ordered_compound_percentage_centered = sweep(ordered_compound_percentage, 2, pca_compound$center, FUN="-")
#Bootstrap PC prediction and reconstruction (61 SNPs -> 7 PCs -> Reconstruct Compounds)
n_bootstrap = 100
ratio_test_data = c(0.1, 0.2, 0.5, 0.8)

bootstrap_results = list()
for (h in 1 : length(ratio_test_data)){
   print(paste0("Bootstrap with test data ratio: ", ratio_test_data[h]))
   bootstrap_results[[h]] = list()
   
  for (i in 1 : n_bootstrap){
    idx_test = sample(1:nrow(pca_compound$x),ceiling(ratio_test_data[h]*nrow(pca_compound$x)),replace=F)
    idx_train = setdiff(1:nrow(pca_compound$x),idx_test)
    
    cur_linear_models = lapply(1:n_pcs, function(j){
      lm(pca_compound$x[idx_train,j] ~ ., as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[j]]),rownames(pca_compound$x)[idx_train],drop=F])))
    })
    predicted_pcs = lapply(1:n_pcs, function(j){
      predict.lm(cur_linear_models[[j]], as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[j]]),rownames(pca_compound$x)[idx_test],drop=F])))
    })
    
    cor_predict_results = sapply(1:n_pcs, function(j){
      cor(pca_compound$x[idx_test,j],predicted_pcs[[j]])
    })
    spearman_predict_results = sapply(1:n_pcs, function(j){
      cor(pca_compound$x[idx_test,j],predicted_pcs[[j]],method="spearman")
    })
    
    bootstrap_results[[h]][[i]] = list(idx_test = idx_test,
                                  predicted_pcs = predicted_pcs,
                                  cor_predict_results = cor_predict_results,
                                  spearman_predict_results = spearman_predict_results)
    print(i)
    
  }
}


#Performance
cor_list = lapply(1:length(bootstrap_results), function(h){
  lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    print(i)
    sapply(1:n_pcs, function(selected_idx){
      predicted_pcs = x$predicted_pcs
      idx_test = x$idx_test
      cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:selected_idx,drop=F] %*% t(pca_compound$rotation[,1:selected_idx,drop=F])
      cur_actual_compound_matrices_centered = ordered_compound_percentage_centered[idx_test, ]
      
      cur_cor_reconstruct_predict = diag(cor(t(cur_reconstruct_result), t(cur_actual_compound_matrices_centered)))
      cur_cor_reconstruct_predict[is.na(cur_cor_reconstruct_predict)] = 0
      return(cur_cor_reconstruct_predict)
    })
  })
})

spearman_list = lapply(1:length(bootstrap_results), function(h){
  lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    print(i)
    sapply(1:n_pcs, function(selected_idx){
      predicted_pcs = x$predicted_pcs
      idx_test = x$idx_test
      cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:selected_idx,drop=F] %*% t(pca_compound$rotation[,1:selected_idx,drop=F])
      cur_actual_compound_matrices_centered = ordered_compound_percentage_centered[idx_test, ]
      
      cur_cor_reconstruct_predict = diag(cor(t(cur_reconstruct_result), t(cur_actual_compound_matrices_centered), method="spearman"))
      cur_cor_reconstruct_predict[is.na(cur_cor_reconstruct_predict)] = 0
      return(cur_cor_reconstruct_predict)
    })
  })
})

  

#LOO Region
loo_region_results = list()
for (i in 1 : length(regions)){
  idx_test = sapply(product_data$`DMD ID`[which(product_data$`Region Code` == i)], function(x){
    which(rownames(pca_compound$x) == x)
  })
  idx_train = setdiff(1:nrow(pca_compound$x),idx_test)
  
  cur_linear_models = lapply(1:n_pcs, function(j){
    lm(pca_compound$x[idx_train,j] ~ ., as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[j]]),rownames(pca_compound$x)[idx_train],drop=F])))
  })
  predicted_pcs = lapply(1:n_pcs, function(j){
    predict.lm(cur_linear_models[[j]], as.data.frame(t(ad_matrix_snp_medium[names(selected_snps_idx[[j]]),rownames(pca_compound$x)[idx_test],drop=F])))
  })
  
  cor_predict_results = sapply(1:n_pcs, function(j){
    cor(pca_compound$x[idx_test,j],predicted_pcs[[j]])
  })
  spearman_predict_results = sapply(1:n_pcs, function(j){
    cor(pca_compound$x[idx_test,j],predicted_pcs[[j]],method="spearman")
  })
  
  loo_region_results[[i]] = list(idx_test = idx_test,
                                predicted_pcs = predicted_pcs,
                                cor_predict_results = cor_predict_results,
                                spearman_predict_results = spearman_predict_results)
  print(i)
  
}

cor_list_loo = lapply(1:length(loo_region_results), function(i){
  x = loo_region_results[[i]]
  print(i)
  sapply(1:n_pcs, function(selected_idx){
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:selected_idx,drop=F] %*% t(pca_compound$rotation[,1:selected_idx,drop=F])
    cur_actual_compound_matrices_centered = sweep(ordered_compound_percentage[idx_test,], 2, pca_compound$center, FUN="-")
    
    cur_cor_reconstruct_predict = diag(cor(t(cur_reconstruct_result), t(cur_actual_compound_matrices_centered)))
    cur_cor_reconstruct_predict[is.na(cur_cor_reconstruct_predict)] = 0
    return(cur_cor_reconstruct_predict)
  })
})

spearman_list_loo = lapply(1:length(loo_region_results), function(i){
  x = loo_region_results[[i]]
  print(i)
  sapply(1:n_pcs, function(selected_idx){
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:selected_idx,drop=F] %*% t(pca_compound$rotation[,1:selected_idx,drop=F])
    cur_actual_compound_matrices_centered = sweep(ordered_compound_percentage[idx_test,], 2, pca_compound$center, FUN="-")
    
    cur_cor_reconstruct_predict = diag(cor(t(cur_reconstruct_result), t(cur_actual_compound_matrices_centered), method="spearman"))
    cur_cor_reconstruct_predict[is.na(cur_cor_reconstruct_predict)] = 0
    return(cur_cor_reconstruct_predict)
  })
})


loo_region_reconstructed_centered = lapply(1:length(loo_region_results), function(i){
  x = loo_region_results[[i]]
  predicted_pcs = x$predicted_pcs
  idx_test = x$idx_test
  cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:n_pcs,drop=F] %*% t(pca_compound$rotation[,1:n_pcs,drop=F])
  cur_actual_compound_matrices_centered = sweep(ordered_compound_percentage[idx_test,], 2, pca_compound$center, FUN="-")
  return(list(reconstructed = cur_reconstruct_result, actual = cur_actual_compound_matrices_centered))
})

library(patchwork)
#Scatter (centered, reconstruct vs actual)


df_scatter_reconstructed_centered = data.frame(
  x = unlist(lapply(loo_region_reconstructed_centered, function(x){as.numeric(as.matrix(x$actual))})),
  y = unlist(lapply(loo_region_reconstructed_centered, function(x){as.numeric(as.matrix(x$reconstructed))})),
  Region = factor(rep(regions, sapply(loo_region_reconstructed_centered, function(x){nrow(x$actual)})), levels=regions)
)

p_scatter_reconstructed_centered = ggplot(df_scatter_reconstructed_centered, aes(x=x, y=y, color=Region)) +
  geom_point_rast(size=1.5, alpha=0.8) +
  theme_classic() +
  labs(title="", x="Actual Compound Profile (Centered)", y="Reconstructed Compound Profile (Centered)") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text = element_text(size=8)) + 
  scale_color_manual(values = brewer.pal(n = 10, name = "Set3")) + 
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

#p_scatter_reconstructed_centered
cor_list_loo_pcs = do.call(rbind, cor_list_loo)
spearman_list_loo_pcs = do.call(rbind, spearman_list_loo)
df_curve_pcs = data.frame(
  y = c(colMeans(cor_list_loo_pcs), colMeans(spearman_list_loo_pcs)),
  ymax = c(colMeans(cor_list_loo_pcs) + 1.96*apply(cor_list_loo_pcs, 2, function(x){sd(x)/sqrt(length(x))}), 
           colMeans(spearman_list_loo_pcs) + 1.96*apply(spearman_list_loo_pcs, 2, function(x){sd(x)/sqrt(length(x))})),
  ymin = c(colMeans(cor_list_loo_pcs) - 1.96*apply(cor_list_loo_pcs, 2, function(x){sd(x)/sqrt(length(x))}), 
           colMeans(spearman_list_loo_pcs) - 1.96*apply(spearman_list_loo_pcs, 2, function(x){sd(x)/sqrt(length(x))})),
  x = rep(1:n_pcs, 2),
  Metric = factor(c(rep("PCC", n_pcs), rep("Spearman", n_pcs)), levels=c("PCC", "Spearman")),
  x2 = sapply(selected_snps_idx, length)[rep(1:n_pcs, 2)]
)
p_curve_pcs = ggplot(df_curve_pcs, aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1.5) + 
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Number of PCs Used for Reconstruction", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text = element_text(size=8)) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = n_pcs),
                     sec.axis = sec_axis(transform = ~.,
                                         breaks = 1:n_pcs,
                                         labels = cumsum(sapply(selected_snps_idx, length)),
                                         name = "Number of SNPs Used for Predicting PCs")) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
#p_curve_pcs

df_boxplot_pcc_region = data.frame(
  Region = factor(sapply(regions, function(i){rep(i,6)}), levels=rev(regions)), 
  PCC = cor_list_loo_pcs[,7]
)
p_boxplot_pcc_region = ggplot(df_boxplot_pcc_region, aes(y=Region, x=PCC, fill=Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0.2, size=1, alpha=0.5) +
  theme_classic() +
  labs(title="", y="Region", x="PCC (Leave One Region Out)") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title.x = element_text(size=8), axis.title.y = element_blank(), axis.text = element_text(size=8)) + 
  scale_fill_manual(values = rev(brewer.pal(n = 10, name = "Set3"))) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))  
p_boxplot_pcc_region

p_final = p_scatter_reconstructed_centered + p_curve_pcs + p_boxplot_pcc_region 
svg(paste0(FIGURE_PATH,"\\SuppFig27_OverallReconstructPerformance.svg"), width = 8, height = 3)
p_final
dev.off()



#Supp data (Curve plot)
actual_loo_region = do.call(rbind, lapply(loo_region_reconstructed_centered, function(x){x$actual}))
reconstructed_loo_region = do.call(rbind, lapply(loo_region_reconstructed_centered, function(x){x$reconstructed}))

pcc_topN = sapply(10:ncol(ordered_compound_percentage), function(j){
  if (j %% 100 == 0) print(j)
  
  cur_reconstructed = reconstructed_loo_region[,1:j]
  cur_actual = actual_loo_region[,1:j]
  return(diag(cor(t(cur_reconstructed), t(cur_actual))))
  
})
colnames(pcc_topN) = 10:ncol(ordered_compound_percentage)

spearman_topN = sapply(10:ncol(ordered_compound_percentage), function(j){
  if (j %% 100 == 0) print(j)
  
  cur_reconstructed = reconstructed_loo_region[,1:j]
  cur_actual = actual_loo_region[,1:j]
  return(diag(cor(t(cur_reconstructed), t(cur_actual), method="spearman")))
})
colnames(spearman_topN) = 10:ncol(ordered_compound_percentage)


pcc_topN_bootstrap = lapply(1:length(bootstrap_results), function(h){
  cur_reconstruct_result = lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:n_pcs,drop=F] %*% t(pca_compound$rotation[,1:n_pcs,drop=F])
    return(cur_reconstruct_result)
  })
  cur_reconstruct_result = do.call(rbind, cur_reconstruct_result)
  
  cur_actual_result = lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_actual_result = ordered_compound_percentage_centered[idx_test,]
    return(cur_actual_result)
  })
  cur_actual_result = do.call(rbind, cur_actual_result)
  
  sapply(c(10,20,50,100,500, 1000, 2000, 3000, 4000, 5000, 5310), function(j){
    if (j %% 100 == 0) print(j)
    cur_reconstruct = cur_reconstruct_result[,1:j]
    cur_actual = cur_actual_result[,1:j]
    return(diag(cor(t(cur_reconstruct), t(cur_actual))))
  })
  
})

spearman_topN_bootstrap = lapply(1:length(bootstrap_results), function(h){
  cur_reconstruct_result = lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_reconstruct_result = do.call(cbind,predicted_pcs)[,1:n_pcs,drop=F] %*% t(pca_compound$rotation[,1:n_pcs,drop=F])
    return(cur_reconstruct_result)
  })
  cur_reconstruct_result = do.call(rbind, cur_reconstruct_result)
  
  cur_actual_result = lapply(1:length(bootstrap_results[[h]]), function(i){
    x = bootstrap_results[[h]][[i]]
    predicted_pcs = x$predicted_pcs
    idx_test = x$idx_test
    cur_actual_result = ordered_compound_percentage_centered[idx_test,]
    return(cur_actual_result)
  })
  cur_actual_result = do.call(rbind, cur_actual_result)
  
  sapply(c(10,20,50,100,500, 1000, 2000, 3000, 4000, 5000, 5310), function(j){
    if (j %% 100 == 0) print(j)
    cur_reconstruct = cur_reconstruct_result[,1:j]
    cur_actual = cur_actual_result[,1:j]
    return(diag(cor(t(cur_reconstruct), t(cur_actual), method="spearman")))
  })
  
})



selected_compound_num = c(10,20,50,100,500, 1000, 2000, 3000, 4000, 5000, 5310)
pcc_topN_selected = pcc_topN[,as.character(selected_compound_num)]
spearman_topN_selected = spearman_topN[,as.character(selected_compound_num)]
df_plot_curve = data.frame(x = c(selected_compound_num, selected_compound_num), 
                           y = c(colMeans(pcc_topN_selected), colMeans(spearman_topN_selected)),
                           ymin = c(colMeans(pcc_topN_selected) - 1.96*apply(pcc_topN_selected, 2, function(x){sd(x)/sqrt(length(x))}), 
                                    colMeans(spearman_topN_selected) - 1.96*apply(spearman_topN_selected, 2, function(x){sd(x)/sqrt(length(x))})),
                           ymax = c(colMeans(pcc_topN_selected) + 1.96*apply(pcc_topN_selected, 2, function(x){sd(x)/sqrt(length(x))}), 
                                    colMeans(spearman_topN_selected) + 1.96*apply(spearman_topN_selected, 2, function(x){sd(x)/sqrt(length(x))})),
                           Metric = factor(c(rep("PCC", length(selected_compound_num)), rep("Spearman", length(selected_compound_num))), levels=c("PCC", "Spearman")))

df_plot_curve_bootstrap = do.call(rbind, lapply(1:length(pcc_topN_bootstrap), function(h){
  data.frame(x = selected_compound_num,
             y = colMeans(pcc_topN_bootstrap[[h]]),
             ymin = colMeans(pcc_topN_bootstrap[[h]]) - 1.96*apply(pcc_topN_bootstrap[[h]], 2, function(x){sd(x)/sqrt(length(x))}), 
             ymax = colMeans(pcc_topN_bootstrap[[h]]) + 1.96*apply(pcc_topN_bootstrap[[h]], 2, function(x){sd(x)/sqrt(length(x))}),
             Metric = "PCC",
             TestDataRatio = paste0("Bootstrap (", ratio_test_data[h]*100, "% Test Data)"))
}))
df_plot_curve_bootstrap = rbind(df_plot_curve_bootstrap, do.call(rbind, lapply(1:length(spearman_topN_bootstrap), function(h){
  data.frame(x = selected_compound_num,
             y = colMeans(spearman_topN_bootstrap[[h]]),
             ymin = colMeans(spearman_topN_bootstrap[[h]]) - 1.96*apply(spearman_topN_bootstrap[[h]], 2, function(x){sd(x)/sqrt(length(x))}), 
             ymax = colMeans(spearman_topN_bootstrap[[h]]) + 1.96*apply(spearman_topN_bootstrap[[h]], 2, function(x){sd(x)/sqrt(length(x))}),
             Metric = "Spearman",
             TestDataRatio = paste0("Bootstrap (", ratio_test_data[h]*100, "% Test Data)"))
})))

df_plot_curve$TestDataRatio = "Leave One Region Out"
final_df_plot_curve = rbind(df_plot_curve, df_plot_curve_bootstrap)
final_df_plot_curve$TestDataRatio = factor(final_df_plot_curve$TestDataRatio, 
                                           levels=c("Leave One Region Out", paste0("Bootstrap (", ratio_test_data*100, "% Test Data)")))


p_curve_loo = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Leave One Region Out"),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(100,500, 1000, 2000, 3000, 4000, 5000)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
p_curve_loo

p_curve_01 = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (10% Test Data)"),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(100,500, 1000, 2000, 3000, 4000, 5000)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_01

p_curve_02 = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (20% Test Data)"),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(100,500, 1000, 2000, 3000, 4000, 5000)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_02

p_curve_05 = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (50% Test Data)"),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text = element_text(size=8)) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(100,500, 1000, 2000, 3000, 4000, 5000)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_05

p_merge = p_curve_loo + p_curve_01 + p_curve_02 + p_curve_05 + plot_layout(ncol=1, guides="collect", axis_titles="collect") & theme(legend.position = "none")



p_curve_looB = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Leave One Region Out" & final_df_plot_curve$x <= 100) ,], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(10, 20, 50, 100)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))

p_curve_01B = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (10% Test Data)" & final_df_plot_curve$x <= 100),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(10, 20, 50, 100)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_01

p_curve_02B = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (20% Test Data)" & final_df_plot_curve$x <= 100),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text.y = element_text(size=8), axis.text.x = element_blank()) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(10, 20, 50, 100)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_02

p_curve_05B = ggplot(final_df_plot_curve[which(final_df_plot_curve$TestDataRatio == "Bootstrap (50% Test Data)" & final_df_plot_curve$x <= 100),], aes(x=x, y=y, color=Metric)) +
  geom_line(linewidth=1) +
  geom_point(size=1.5) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=Metric), alpha=0.2, color=NA) +
  theme_classic() +
  labs(title="", x="Top N Compounds", y="Correlation") +
  theme(plot.title = element_text(hjust = 0.5, size=10), legend.text = element_text(size=8), legend.spacing.y=unit(-0.3,"cm"), legend.position = "none",
        axis.title = element_text(size=8), axis.text = element_text(size=8)) + 
  scale_color_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) + 
  scale_fill_manual(values=c("PCC" = "steelblue", "Spearman" = "#E41A1C")) +
  scale_x_continuous(breaks = c(10, 20, 50, 100)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,0.8))
#p_curve_05

svg(paste0(FIGURE_PATH,"/SuppFig28_ReconstructPerformance_Bootstrap.svg"), width = 8, height = 8)
p_merge = (p_curve_loo + p_curve_looB) / (p_curve_01 + p_curve_01B) / (p_curve_02 + p_curve_02B) / (p_curve_05 + p_curve_05B) + 
  plot_layout(guides="collect", axis_titles="collect") & 
  theme(plot.margin = margin(2, 2, 2, 2, "pt"), legend.position = "bottom")
print(p_merge)
dev.off()


