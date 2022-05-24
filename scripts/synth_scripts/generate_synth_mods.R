prepare_sim_data <- function(){
  
  #Prepare data
  #Geneinfo
  gene_info_p <- paste0(getwd(),'/prep_simulated/promote_v1.gene_info.txt')
  gi <- read.delim(gene_info_p)[,c('uniq_isos','regulator')]
  rownames(gi) <- gi[,1]
  
  #exp
  exp_p <- paste0(getwd(),'/prep_simulated/mrm_sm_3500.txt')
  exp <- read.delim(exp_p)
  
  #linkeroutput
  
  #clinic file
  clinic_p <- paste0(getwd(),'/prep_simulated/promote_clinical_ours_fixed.txt')
  clinic <- read.delim(clinic_p)
  rownames(clinic) <- clinic[,1]
  
  #Generate regulatordata
  regs <- gi[which(gi$regulator==1),1]
  regData <- exp[regs,]
  
  return(list(regs=regData,gene_p=gene_info_p,exp_p=exp_p,clinic_p=clinic_p,clinic=clinic))
}

generate_sim_data<-function(regulatorData,clinic,NrModules=10,
                            lambda=0.15,dorew=FALSE,specificrew=c()){
  
  
  #Get the regulator names
  rnames <- rownames(regulatorData)
  
  #Convert to matrix
  if (inherits(regulatorData,'data.frame')){
    regulatorData <- as.matrix(regulatorData)
  }
  
  #Get dimensions
  NrSamples<-ncol(regulatorData)
  NrRegulators<-nrow(regulatorData)
  
  maxNrGenes <- 1e6
  
  ExpectedNrGenesPerModule<-200
  ExpectedNrRegulatorsPerModule<-5
  ExpectedRewiredModules <- 3
  
  if (dorew){
    #Choose a % of rewired modules to introduce the difference of noise
    if (!length(specificrew)){
      rewmods <- rbinom(NrModules,1,ExpectedRewiredModules/NrModules)
    }else{
      rewmods <- rep(0,NrModules)
      rewmods[specificrew] <- 1
    }
    message('Rewired modules ',sum(rewmods))
  }
  
  #Initialize the list of modules
  ModulesData<-list()
  #Initialize the regulators 
  modulesRegulators<-matrix(data=NA,nrow=NrRegulators, ncol=NrModules)
  rownames(modulesRegulators) <- rownames(regulatorData)
  #Initialize the number of target genes we have "generated"
  tot_target_genes <- 0
  
  grp <- as.numeric(sapply(colnames(regulatorData),function(x) clinic[x,2]))
  #Define NR (0) as we are going to increase the noise in that part
  clinic_NR_pos <- which(grp==0)
  #Define samples sorting clinic file by phenotype (First Responders (1), then NonResponders (0))
  clinic_names <- clinic$Sample.ID[order(clinic$Class,decreasing = TRUE)]
  
  for(i in 1:NrModules){
    
    #Initialize empty vector
    moduleRegulators<-numeric(NrRegulators)
    
    #Sample (at least 1 regulator)
    while(sum(moduleRegulators)<=1){
      moduleRegulators<-rbinom(NrRegulators, 1, ExpectedNrRegulatorsPerModule/NrRegulators)
    }
    
    #Generate the module mean as a linear combinations of randomly sampled drivers
    regulators<-which(moduleRegulators==1)
    #Calculate the beta of each regulator
    betas<-rnorm(length(regulators), 0,1)
    #Store it
    moduleRegulators[regulators]<-betas
    modulesRegulators[,i]<-moduleRegulators
    
    #Compute the linear comb
    module_mean<-as.vector(colSums(moduleRegulators*regulatorData))
    #Scale it
    module_mean<-as.vector(scale(module_mean))
    
    #message('The module mean is: ',mean(module_mean),'\n')
    
    #Generate the number of target genes in the module
    NrGenesInModule<-rbinom(1, maxNrGenes, ExpectedNrGenesPerModule/maxNrGenes)
    
    #Add the noise
    modmean_l <- length(module_mean)
    
    
    for (lambda_e in lambda){
      
      #method 1, rnorm for each gene and sample + modify the mean by sample
      #As the variance is very low and they come from the same rnorm, the will have high covariance
      geneModule <- t(module_mean+t(matrix(rnorm(NrGenesInModule*modmean_l,mean=0,sd=lambda_e), 
                                           NrGenesInModule,modmean_l)))
      #scale the genes
      geneModule <- t(scale(t(geneModule)))
      
      if (rewmods[i] & dorew){
        
        #For the rewired, we do the whitening on the NR samples
        old_NR_module <- geneModule[,clinic_NR_pos]
        #scale it 
        old_NR_module <- t(scale(t(old_NR_module)))
        #calculate cov
        cov_old <- t(cov(t(old_NR_module)))
        #do eigen
        eigen_cov_old_D <- diag(abs(eigen(cov_old)$values))
        eigen_cov_old_tU <- t(eigen(cov_old)$vectors)
        
        #apply the whitening by y = inv(sqrt(D)) * t(U) * x
        new_NR_module <- solve(sqrt(eigen_cov_old_D))%*%eigen_cov_old_tU%*%old_NR_module
        
        geneModule[,clinic_NR_pos] <- new_NR_module
        
        #scale the genes
        geneModule <- t(scale(t(geneModule)))
        
        
      }
      
      rownames(geneModule) <- paste0('g',tot_target_genes+seq(nrow(geneModule)))
      
      #Update the number of target genes
      tot_target_genes <- tot_target_genes + NrGenesInModule
      
      #Add the drivers
      fullModule <- rbind(regulatorData[regulators,],geneModule)
      
      ModulesData[[paste0(c('mod',i,lambda_e),collapse='_')]] <- fullModule
      
    }
    
    #Check if maximum number of target_genes has been reached
    #if(tot_target_genes>maxNrGenes) break
    
  }
  
  #Create our final list of synthetic datasets
  synth_list <- list()
  
  
  for (lambda_e in lambda){
    
    grep_lambda <- paste0('mod_','[0-9]+_',lambda_e,'$')
    lambda_modules <- grep(grep_lambda,names(ModulesData))
    
    lambda_e_ModulesData <- ModulesData[lambda_modules]
    message(length(lambda_e_ModulesData),' modules generated with noise ',lambda_e)
    
    #From list of modules to exp matrix
    drivers_mod <- sapply(lambda_e_ModulesData,
                          function(x) intersect(rownames(regulatorData),rownames(x)))
    synth_drivers <- unique(unlist(drivers_mod))
    
    synth_targs <- lapply(lambda_e_ModulesData,function(x) x[setdiff(rownames(x),synth_drivers),])
    
    #Concat rows of synth_targs
    synth_targs_matrix <- do.call(rbind,synth_targs)
    
    #Add the regs
    expmat <- as.matrix(rbind(regulatorData[synth_drivers,], synth_targs_matrix))
    
    #Set colnames 
    colnames(expmat) <- clinic_names
    
    synth_list[[paste0(c('noise',lambda_e),collapse='_')]] <- list(modules=lambda_e_ModulesData,expmat=expmat,
                                                                   regs=synth_drivers,targs=rownames(synth_targs_matrix),
                                                                   drivers_mod=drivers_mod,regprograms=t(modulesRegulators),
                                                                   rewmods = rewmods)
  }
  
  return(synth_list)
  
}

generate_sim_data_samples<-function(regulatorData_orig,clinic_orig,NrModules=10,
                                    lambda=0.15,dorew=FALSE,specificrew=c(),nrsamples_v=FALSE){
  
  #Create our final list of synthetic datasets
  synth_list <- list()
  
  for (nrsamples in nrsamples_v){
    
    message('Nr of samples ',nrsamples)
    
    #divide clinic in phenotype
    Rsamples <- which(clinic_orig$Class==1)
    NRsamples <- which(clinic_orig$Class==0)
    
    #Calculate min to maintain proportions
    prop_min <- min(length(NRsamples),length(NRsamples))
    
    if ((nrsamples/2)>prop_min){
      if (length(Rsamples)>length(NRsamples)){
        ##Sample maintaining proportions
        sampled_samples <- c(sample(Rsamples,prop_min+nrsamples-(prop_min*2)),sample(NRsamples,prop_min))
      }else{
        ##Sample maintaining proportions
        sampled_samples <- c(sample(Rsamples,prop_min),sample(NRsamples,nrsamples-(prop_min*2)))
      }
      
    }else{
      ##Sample maintaining proportions
      sampled_samples <- c(sample(Rsamples,nrsamples/2),sample(NRsamples,nrsamples/2))
    }
    
    
    
    #Filter clinic file
    clinic <- clinic_orig[sampled_samples,]
    regulatorData <- regulatorData_orig[,clinic$Sample.ID]
    
    ##Finish sampling
    
    #Get the regulator names
    rnames <- rownames(regulatorData)
    
    #Convert to matrix
    if (inherits(regulatorData,'data.frame')){
      regulatorData <- as.matrix(regulatorData)
    }
    
    #Get dimensions
    NrSamples<-ncol(regulatorData)
    NrRegulators<-nrow(regulatorData)
    
    maxNrGenes <- 1e6
    
    ExpectedNrGenesPerModule<-200
    ExpectedNrRegulatorsPerModule<-5
    ExpectedRewiredModules <- 3
    
    if (dorew){
      #Choose a % of rewired modules to introduce the difference of noise
      if (!length(specificrew)){
        rewmods <- rbinom(NrModules,1,ExpectedRewiredModules/NrModules)
      }else{
        rewmods <- rep(0,NrModules)
        rewmods[specificrew] <- 1
      }
      message('Rewired modules ',sum(rewmods))
    }
    
    #Initialize the list of modules
    ModulesData<-list()
    #Initialize the regulators 
    modulesRegulators<-matrix(data=NA,nrow=NrRegulators, ncol=NrModules)
    rownames(modulesRegulators) <- rownames(regulatorData)
    #Initialize the number of target genes we have "generated"
    tot_target_genes <- 0
    
    grp <- as.numeric(sapply(colnames(regulatorData),function(x) clinic[x,2]))
    #Define NR (0) as we are going to increase the noise in that part
    clinic_NR_pos <- which(grp==0)
    #Define samples sorting clinic file by phenotype (First Responders (1), then NonResponders (0))
    clinic_names <- clinic$Sample.ID[order(clinic$Class,decreasing = TRUE)]
    
    for(i in 1:NrModules){
      
      #Initialize empty vector
      moduleRegulators<-numeric(NrRegulators)
      
      #Sample (at least 1 regulator)
      while(sum(moduleRegulators)<=1){
        moduleRegulators<-rbinom(NrRegulators, 1, ExpectedNrRegulatorsPerModule/NrRegulators)
      }
      
      #Generate the module mean as a linear combinations of randomly sampled drivers
      regulators<-which(moduleRegulators==1)
      #Calculate the beta of each regulator
      betas<-rnorm(length(regulators), 0,1)
      #Store it
      moduleRegulators[regulators]<-betas
      modulesRegulators[,i]<-moduleRegulators
      
      #Compute the linear comb
      module_mean<-as.vector(colSums(moduleRegulators*regulatorData))
      #Scale it
      module_mean<-as.vector(scale(module_mean))
      
      #message('The module mean is: ',mean(module_mean),'\n')
      
      #Generate the number of target genes in the module
      NrGenesInModule<-rbinom(1, maxNrGenes, ExpectedNrGenesPerModule/maxNrGenes)
      
      #Add the noise
      modmean_l <- length(module_mean)
      
      #method 1, rnorm for each gene and sample + modify the mean by sample
      #As the variance is very low and they come from the same rnorm, the will have high covariance
      geneModule <- t(module_mean+t(matrix(rnorm(NrGenesInModule*modmean_l,mean=0,sd=lambda), 
                                           NrGenesInModule,modmean_l)))
      #scale the genes
      geneModule <- t(scale(t(geneModule)))
      
      if (rewmods[i] & dorew){
        
        #For the rewired, we do the whitening on the NR samples
        old_NR_module <- geneModule[,clinic_NR_pos]
        #scale it 
        old_NR_module <- t(scale(t(old_NR_module)))
        #calculate cov
        cov_old <- t(cov(t(old_NR_module)))
        #do eigen
        eigen_cov_old_D <- diag(abs(eigen(cov_old)$values))
        eigen_cov_old_tU <- t(eigen(cov_old)$vectors)
        
        #apply the whitening by y = inv(sqrt(D)) * t(U) * x
        new_NR_module <- solve(sqrt(eigen_cov_old_D))%*%eigen_cov_old_tU%*%old_NR_module
        
        geneModule[,clinic_NR_pos] <- new_NR_module
        
        #scale the genes
        geneModule <- t(scale(t(geneModule)))
        
        
      }
      
      rownames(geneModule) <- paste0('g',tot_target_genes+seq(nrow(geneModule)))
      
      #Update the number of target genes
      tot_target_genes <- tot_target_genes + NrGenesInModule
      
      #Add the drivers
      fullModule <- rbind(regulatorData[regulators,],geneModule)
      
      ModulesData[[paste0(c('mod',i,nrsamples),collapse='_')]] <- fullModule
      
      
      
      #Check if maximum number of target_genes has been reached
      #if(tot_target_genes>maxNrGenes) break
      
    }
    
    grep_samples <- paste0('mod_','[0-9]+_',nrsamples,'$')
    samples_modules <- grep(grep_samples,names(ModulesData))
    
    samples_ModulesData <- ModulesData[samples_modules]
    message(length(samples_ModulesData),' modules generated with samples ',nrsamples)
    
    #From list of modules to exp matrix
    drivers_mod <- sapply(samples_ModulesData,
                          function(x) intersect(rownames(regulatorData),rownames(x)))
    synth_drivers <- unique(unlist(drivers_mod))
    
    synth_targs <- lapply(samples_ModulesData,function(x) x[setdiff(rownames(x),synth_drivers),])
    
    #Concat rows of synth_targs
    synth_targs_matrix <- do.call(rbind,synth_targs)
    
    #Add the regs
    expmat <- as.matrix(rbind(regulatorData[synth_drivers,], synth_targs_matrix))
    
    #Set colnames 
    colnames(expmat) <- clinic_names
    
    synth_list[[paste0(c('samples',nrsamples),collapse='_')]] <- list(modules=samples_ModulesData,expmat=expmat,
                                                                      regs=synth_drivers,targs=rownames(synth_targs_matrix),
                                                                      drivers_mod=drivers_mod,regprograms=t(modulesRegulators),
                                                                      rewmods = rewmods)
  }
  
  return(synth_list)
  
}


#Prepare the data
prep_sim <- prepare_sim_data()

#Generating synths
lambda_noise <- c(0.01,seq(0.4,5,2/5),10,25,50,100)

for (lambda in lambda_noise){
  #Call synthethics
  synth_modules <- generate_sim_data(prep_sim$regs,prep_sim$clinic,NrModules=10,
                                     lambda=lambda,dorew=TRUE,specificrew=c(1,4,7))
  
  saveRDS(synth_modules,paste0(getwd(),'/input_data/from_0.01_to_100_noise_rewired/synth_modules_rewired_noise_',lambda,'.rds'))
}

#Now lets generate synthetic data but varying the number of samples
for (num_samples in c(10,16,28,38,46)){
  df <- generate_sim_data_samples(regulatorData_orig=prep_sim$regs,clinic_orig=prep_sim$clinic,NrModules=10,
                            lambda=1.2,dorew=TRUE,specificrew = c(1,4,7),nrsamples_v=num_samples)
  saveRDS(df,file=paste0(getwd(),'/input_data/from_18_to_46_samples_rewired/synth_module_',num_samples,'_samples.rds'))
}
