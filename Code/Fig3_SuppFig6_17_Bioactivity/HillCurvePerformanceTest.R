#Definition of the path:
#1. Working directory (PATH)
#2. Data directory (DATA_PATH)
#3. Processed data directory (PROCESSED_DATA_PATH)
if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
source(paste0(this.path::this.dir(),"\\..\\FilePath.R"))
source(paste0(CODE_PATH,"\\Analysis_Utils.r"))


if (!file.exists(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"))){ 
  stop("A representative Hill curve evaluated based on the small molecules dataset is needed. Please run the script 'SmallMoleculeBioactivityEfficacyEstimation.r' first to extract the Hill curves (MilkSmallMoleculeBioactivitySummary.csv).")
}

library(stringr)

curve_table = read.csv(paste0(PROCESSED_DATA_PATH,"\\Fig3_SuppFig6_17_Bioactivity\\MilkSmallMoleculeBioactivitySummary.csv"), header = T, check.names = F)
curve_table = curve_table[which(curve_table[,"Fit_ZeroActivity"] != curve_table[,"Fit_InfiniteActivity"]),]

library(basicdrm)

get_representative_hill_curve = function(curve_table_selected, fold_change = 10^seq(-5,5,0.1)){
  par_table = curve_table_selected[,c("Fit_LogAC50","Fit_HillSlope","Fit_ZeroActivity","Fit_InfiniteActivity")]
  par_table[,1] = 10^(par_table[,1]) #Convert logAC50 to AC50
  colnames(par_table) = c("AC50","HillSlope","ZeroActivity","InfiniteActivity")
  
  efficacy_eval = sapply(fold_change, function(x){
    cur_efficacy = sapply(1:nrow(par_table), function(i){
      hpar = as.numeric(par_table[i,])
      return(basicdrm::evalHillModel(x*hpar[1], hpar))
    })
    return(cur_efficacy)
  })
  colnames(efficacy_eval) = log10(fold_change)
  
  efficacy_eval_converted = apply(efficacy_eval,2,function(x){
    (x - par_table[,"ZeroActivity"])/(par_table[,"InfiniteActivity"] - par_table[,"ZeroActivity"])
  })
  ratio_efficacy_eval = efficacy_eval_converted/0.5
  
  ratio_efficacy_eval_mean = apply(ratio_efficacy_eval,2,mean)
  ratio_efficacy_eval_mean_curve = fitHillModel(
    10^seq(-5,5,0.1),
    as.numeric(ratio_efficacy_eval_mean),
    c(1,0.5,0,2),
    start=c(1,0.5,0,2),
  )
  
  return(ratio_efficacy_eval_mean_curve)
}


#LOO
fold_change = 10^seq(-5,5,0.1)
par_truth = NULL
par_representative_loo = NULL

val_truth = NULL
val_representative_loo = NULL

for (i in 1 : nrow(curve_table)){
  curve_table_loo = curve_table[-i,]
  representative_curve_loo = get_representative_hill_curve(curve_table_loo)
  
  par_representative_loo = rbind(par_representative_loo, as.numeric(representative_curve_loo$coefficients))
  cur_par_truth = as.numeric(curve_table[i,c("Fit_LogAC50","Fit_HillSlope","Fit_ZeroActivity","Fit_InfiniteActivity")])
  cur_par_truth[1] = 1
  cur_par_truth[3] = 0
  cur_par_truth[4] = 2
  
  par_truth = rbind(par_truth, cur_par_truth)
  
  cur_val_truth = basicdrm::evalHillModel(fold_change, cur_par_truth)
  cur_val_representative_loo = basicdrm::evalHillModel(fold_change, as.numeric(representative_curve_loo$coefficients))
  
  val_truth = rbind(val_truth, cur_val_truth * 50) #convert to 0 ~ 100%
  val_representative_loo = rbind(val_representative_loo, cur_val_representative_loo * 50) #convert to 0 ~ 100%
}

colnames(val_truth) = seq(-5,5,0.1)
colnames(val_representative_loo) = seq(-5,5,0.1)



#LOO Graph
idx_shown = 51:81
df_estimation_error = data.frame(
  avg_error = as.vector(colMeans(abs(val_representative_loo[,idx_shown] - val_truth[,idx_shown]))),
  error_upper = as.vector(apply((abs(val_representative_loo[,idx_shown] - val_truth[,idx_shown])),2,function(x){
    mean(x) + sd(x)*2
  })),
  error_lower = as.vector(apply((abs(val_representative_loo[,idx_shown] - val_truth[,idx_shown])),2,function(x){
    mean(x) - sd(x)*2
  })),
  error = as.vector(abs(val_representative_loo[,idx_shown] - val_truth[,idx_shown])),
  fold_change = round(as.numeric(colnames(val_truth)[idx_shown]),1)
)

library(ggplot2)


#Bootstrap
run_bootstrap = function(n_bootstrap, ratio_select, seed){
  val_representative_bootstrap = array(0, dim = c(n_bootstrap, nrow(curve_table) - (ceiling(nrow(curve_table) * ratio_select)), length(fold_change)))
  val_representative_bootstrap_truth = array(0, dim = c(n_bootstrap, nrow(curve_table) - (ceiling(nrow(curve_table) * ratio_select)), length(fold_change)))
  set.seed(seed)
  
  for (i in 1 : n_bootstrap){
    print(i)
    sample_idx = sample(1:nrow(curve_table), ceiling(nrow(curve_table) * ratio_select), replace = F)
    curve_table_train = curve_table[sample_idx,]
    curve_table_test = curve_table[setdiff(1:nrow(curve_table),sample_idx),]
    representative_curve_bootstrap = get_representative_hill_curve(curve_table_train)
    for (j in 1 : nrow(curve_table_test)){
      cur_par_truth = as.numeric(curve_table_test[j,c("Fit_LogAC50","Fit_HillSlope","Fit_ZeroActivity","Fit_InfiniteActivity")])
      cur_par_truth[1] = 1
      cur_par_truth[3] = 0
      cur_par_truth[4] = 2
      cur_val_truth = basicdrm::evalHillModel(fold_change, cur_par_truth)
      val_representative_bootstrap_truth[i,j,] = cur_val_truth * 50 #convert to 0 ~ 100%
      val_representative_bootstrap[i,j,] = basicdrm::evalHillModel(fold_change, as.numeric(representative_curve_bootstrap$coefficients)) * 50
    }
  }
  return(list(
    val_representative_bootstrap = val_representative_bootstrap,
    val_representative_bootstrap_truth = val_representative_bootstrap_truth
  ))
}
run_bootstrap_09 = run_bootstrap(n_bootstrap = 100, ratio_select = 0.9, seed = 123)
run_bootstrap_05 = run_bootstrap(n_bootstrap = 100, ratio_select = 0.5, seed = 123)
run_bootstrap_01 = run_bootstrap(n_bootstrap = 100, ratio_select = 0.1, seed = 123)


idx_shown = 51:81
get_df_estimation_error_bootstrap = function(run_bootstrap, idx_shown){
  val_representative_bootstrap = run_bootstrap$val_representative_bootstrap
  val_representative_bootstrap_truth = run_bootstrap$val_representative_bootstrap_truth
  
  df_estimation_error_bootstrap = data.frame(
    avg_error = rowMeans(sapply(1:dim(val_representative_bootstrap)[1],function(i){
      colMeans(abs(val_representative_bootstrap[i,,idx_shown] - val_representative_bootstrap_truth[i,,idx_shown]))
    })),
    error_upper95 = apply(sapply(1:dim(val_representative_bootstrap)[1],function(i){
      apply((abs(val_representative_bootstrap[i,,idx_shown] - val_representative_bootstrap_truth[i,,idx_shown])),2,function(x){quantile(x,0.95)})
    }),1,mean),
    error_upper = apply(sapply(1:dim(val_representative_bootstrap)[1],function(i){
      apply((abs(val_representative_bootstrap[i,,idx_shown] - val_representative_bootstrap_truth[i,,idx_shown])),2,max)
    }),1,max),
    error_lower = apply(sapply(1:dim(val_representative_bootstrap)[1],function(i){
      apply((abs(val_representative_bootstrap[i,,idx_shown] - val_representative_bootstrap_truth[i,,idx_shown])),2,min)
    }),1,min),
    fold_change = round(seq(-5,5,0.1)[idx_shown],1)
  )
  
  return(df_estimation_error_bootstrap)
}

df_estimation_error_bootstrap_09 = get_df_estimation_error_bootstrap(run_bootstrap_09, idx_shown)
df_estimation_error_bootstrap_05 = get_df_estimation_error_bootstrap(run_bootstrap_05, idx_shown)
df_estimation_error_bootstrap_01 = get_df_estimation_error_bootstrap(run_bootstrap_01, idx_shown)

library(ggplot2)
library(RColorBrewer)
svg(paste0(FIGURE_PATH,"\\SuppFig17b_EstimatedHillCurvePerformance.svg"), width = 4, height = 4)
h = ggplot(data = df_estimation_error_bootstrap_09) + 
  geom_line(aes(x = fold_change, y = avg_error),color="#66C2A5", linewidth=1) +
  geom_point(aes(x = fold_change, y = avg_error),color="#66C2A5") +
  geom_line(aes(x = fold_change, y = error_upper95),linetype="dashed",color="#66C2A5", linewidth=1) +
  #geom_point(aes(x = fold_change, y = error_upper95),color="#66C2A5") +
  #geom_ribbon(data = df_estimation_error_bootstrap_09, aes(x = fold_change, ymin = 0, ymax = error_upper), alpha = 0.5, fill="#66C2A5") +
  #geom_line(aes(x = fold_change, y = error_upper),color="#66C2A5") +
  
  geom_line(data = df_estimation_error_bootstrap_05, aes(x = fold_change, y = avg_error),color="#FC8D62", linewidth=1) +
  geom_point(data = df_estimation_error_bootstrap_05,aes(x = fold_change, y = avg_error),color="#FC8D62") +
  geom_line(data = df_estimation_error_bootstrap_05, aes(x = fold_change, y = error_upper95),linetype="dashed",color="#FC8D62", linewidth=1) +
  #geom_point(data = df_estimation_error_bootstrap_05, aes(x = fold_change, y = error_upper95),color="#FC8D62") +
  #geom_ribbon(data = df_estimation_error_bootstrap_05, aes(x = fold_change, ymin = 0, ymax = error_upper), alpha = 0.5, fill="#FC8D62") +
  #geom_line(data = df_estimation_error_bootstrap_05, aes(x = fold_change, y = error_upper),color="#FC8D62") +
  
  geom_line(data = df_estimation_error_bootstrap_01, aes(x = fold_change, y = avg_error),color="#8DA0CB", linewidth=1) +
  geom_point(data = df_estimation_error_bootstrap_01, aes(x = fold_change, y = avg_error),color="#8DA0CB") +
  geom_line(data = df_estimation_error_bootstrap_01, aes(x = fold_change, y = error_upper95),linetype="dashed",color="#8DA0CB", linewidth=1) +
  #geom_point(data = df_estimation_error_bootstrap_01, aes(x = fold_change, y = error_upper95),color="#8DA0CB") +
  #geom_ribbon(data = df_estimation_error_bootstrap_01, aes(x = fold_change, ymin = 0, ymax = error_upper), alpha = 0.5, fill="#8DA0CB") +
  #geom_line(data = df_estimation_error_bootstrap_01, aes(x = fold_change, y = error_upper),color="#8DA0CB") +
  
  xlab("Abs Log Fold Change from IC50") + ylab("Absolute Error in Efficacy Estimation (%)") +
  theme_classic() + 
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30,5))
print(h)
dev.off()









