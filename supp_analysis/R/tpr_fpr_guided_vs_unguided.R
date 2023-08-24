library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#setwd("~/Documents/GSFA/")
#setwd("/project2/xinhe/yifan/Factor_analysis/simulations/simulation_results/multi_factors_n_markers/realistic_count_beta_0.1_to_0.6/unguided_vs_guided.pi_0.05/single_datasets/")
setwd("/project2/xinhe/yifan/Factor_analysis/simulations/simulation_results/multi_factors_n_markers/realistic_beta_0.1_to_0.6/unguided_vs_guided.pi_0.05/single_datasets/")

fnr.guided <- matrix(0, nrow=300, ncol=6, dimnames = list(paste0("rep",1:300), paste0("guide",1:6)))
fpr.guided <- matrix(0, nrow=300, ncol=6, dimnames = list(paste0("rep",1:300), paste0("guide",1:6)))
fnr.unguided <- matrix(0, nrow=300, ncol=6, dimnames = list(paste0("rep",1:300), paste0("guide",1:6)))
fpr.unguided <- matrix(0, nrow=300, ncol=6, dimnames = list(paste0("rep",1:300), paste0("guide",1:6)))

for (rep in 1:300) {
  cat("working on rep", rep)
  fit_result <- readRDS(paste0("fit_result.rep_",rep,".rds"))
  
  deg.factor <- foreach(i=1:10) %do% which(fit_result$unguided$posterior_means$F_pm[,i] > 0.95)
  
  temp <- as.data.frame(cbind(fit_result$sim_data$G, fit_result$unguided$posterior_means$Z_pm))
  colnames(temp) <- c(paste0("guide",1:6), paste0("factor",1:10))
  beta.pval <- matrix(0, nrow=6, ncol=10)
  for(f in 1:10) {
    beta.pval[,f] <- foreach (i=1:6, .combine = c) %dopar% {
      form.cur <- as.formula(paste0("factor",f, "~guide", i))
      res <- lm(form.cur, data=temp)
      summary(res)$coefficients[2,4]
    }
  }
  
  beta.qval <- matrix(p.adjust(c(beta.pval),"BH"), nrow=6)
  proxies <- foreach(i=1:6) %do% which(beta.qval[i,]<0.05)
  deg.guide <- foreach(i=1:6) %do% unique(unlist(deg.factor[proxies[[i]]]))
  theta.true <- (fit_result$sim_data$U * fit_result$sim_data$F) %*% t(fit_result$params$beta_true)
  fpr.unguided[rep,] <- foreach(i=1:6, .combine=c) %do% {
    sum(deg.guide[[i]] %in% which(theta.true[,i]==0)) / sum(theta.true[,i]==0)
  }
  fnr.unguided[rep,] <- 1 - foreach(i=1:6, .combine=c) %do% {
    sum(deg.guide[[i]] %in% which(abs(theta.true[,i])>0)) / sum(abs(theta.true[,i])>0)
  }
  fpr.guided[rep,] <- foreach(i=1:6, .combine=c) %do% {
    sum(which(fit_result$GSFA$lfsr[,i]<.05) %in% which(theta.true[,i]==0)) / sum(theta.true[,i]==0)
  }
  fnr.guided[rep,] <- 1 - foreach(i=1:6, .combine=c) %do% {
    sum(which(fit_result$GSFA$lfsr[,i]<.05) %in% which(abs(theta.true[,i])>0)) / sum(abs(theta.true[,i])>0)
  }
}
#save(fpr.unguided,fpr.guided, fnr.guided, fnr.unguided, file="fit_count_fpr_fnr.Rdata")
save(fpr.unguided,fpr.guided, fnr.guided, fnr.unguided, file="/project2/xinhe/lifanl/GSFA_sim/fit_normal_fpr_fnr.Rdata")


