#Check rewiring
check_rewiring_synth <- function(synth_modules,clinic){
  
  #first select rewired and not rewired
  rew_pos <- which(synth_modules$rewmods==1)[1]
  
  #we are going to check that the score between a rewired and not rewired is different
  rewired_lognorm <- synth_modules$modules[[rew_pos]]
  
  grp <- as.numeric(sapply(colnames(rewired_lognorm),function(x) clinic[x,2]))
  #compute cov matrices
  NR_pos <- which(grp==0)
  R_pos <- which(grp==1)
  
  #of rewired
  rew_cov_R <- t(cor(t(rewired_lognorm[,R_pos])))
  rew_cov_NR <- t(cor(t(rewired_lognorm[,NR_pos])))
  
  #calculate score (euclidean distance)
  edist_rew <- sum((rew_cov_R-rew_cov_NR)^2)
  
  edist_random <- c()
  #now do random sampling and compute the euclidean distance
  
  for (i in seq(1e3)){
    
    random_samples <- sample(seq(nrow(clinic)),nrow(clinic))
    random_int <- sample(seq(5,nrow(clinic)-5),1)
    #of rewired
    rew_cov_R <- t(cor(t(rewired_lognorm[,random_samples[1:random_int]])))
    rew_cov_NR <- t(cor(t(rewired_lognorm[,random_samples[random_int:nrow(clinic)]])))
    
    edist_random <- c(edist_random,sum((rew_cov_R-rew_cov_NR)^2))
    
  }
  edist_random <- edist_random[!is.na(edist_random)]
  
  message('Score with right choice ',edist_rew)
  message('95% quant score with random choice ', quantile(edist_random,probs=0.95))
  
}
check_rewiring_linker_homemade <- function(lognorm,clinic,linkerout,mod){
  
  #Get drivers and targets of the rewired mod
  regs <- linkerout$modules[[1]][[mod]]$regulators
  targs <- linkerout$modules[[1]][[mod]]$target_genes
  
  #Retrieve the rewired matrix
  rewired_lognorm <- lognorm[c(regs,targs),]
  
  #compute cov matrices
  NR_pos <- which(clinic$Class==0)
  R_pos <- which(clinic$Class==1)
  
  #of rewired
  rew_cov_R <- t(cor(t(rewired_lognorm[,R_pos])))
  rew_cov_NR <- t(cor(t(rewired_lognorm[,NR_pos])))
  
  #calculate score (euclidean distance)
  edist_rew <- sum((rew_cov_R-rew_cov_NR)^2)
  
  message('Score is: ',edist_rew)
  
  edist_random <- c()
  #now do random sampling and compute the euclidean distance
  
  for (i in seq(1e3)){
    
    random_samples <- sample(seq(nrow(clinic)),nrow(clinic))
    random_int <- sample(seq(5,nrow(clinic)-5),1)
    #of rewired
    rew_cov_R <- t(cor(t(rewired_lognorm[,random_samples[1:random_int]])))
    rew_cov_NR <- t(cor(t(rewired_lognorm[,random_samples[random_int:nrow(clinic)]])))
    
    edist_random <- c(edist_random,sum((rew_cov_R-rew_cov_NR)^2))
    
  }
  edist_random <- edist_random[!is.na(edist_random)]
  
  message('Score with right choice ',edist_rew)
  message('95% quant score with random choice ', quantile(edist_random,probs=0.95))
  
}
check_rewiring_linker_dave <- function(x, clinic, perm = 500){
  
  grp <- as.numeric(sapply(rownames(x),function(x) clinic[x,2]))+1
  
  p = ncol(x)
  numgrp = length(unique(grp))
  
  ## keep only rows with variance within groups
  bools = rep(TRUE, p)
  for (g in unique(grp)) {
    bools = bools & (matrixStats::colSds(x[grp == g, ]) != 0)
  }
  if (sum(bools) < p) {
    methods::show(paste(collapse = " ", c("...Rewiring Test - Dropping 0 Variance Genes:", paste(collapse = ",", names(which(bools == FALSE))))))
    # show(c(p,sum(bools),which(bools==FALSE)))
    x = x[, bools]
    p = ncol(x)
  }
  
  ## test stat
  T = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
    as.vector((stats::cor(x[grp == pair[1], seq_len(p)]) - stats::cor(x[grp == pair[2], seq_len(p)]))^2)
  }))
  
  ## p
  T_star = rep(NA, perm)
  for (j in seq_len(perm)) {
    grp_perm = sample(grp)
    T_star[j] = sum(apply(utils::combn(numgrp, 2), 2, function(pair) {
      as.vector((stats::cor(x[grp_perm == pair[1], ]) - stats::cor(x[grp_perm == pair[2], ]))^2)
    }))
  }
  # show(c(T, T_star))
  return(mean(c(T, T_star) >= T, na.rm = TRUE))
  
}
