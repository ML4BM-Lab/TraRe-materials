difmodels_difnoise <- function(synth_modules, methodlink='LASSOmin',FDR=0.05){
  
  message('Working with.. ',methodlink)
  
  for (synth_modules_name in names(synth_modules)){
    
    lambda_noise <- as.numeric(substr(synth_modules_name,7,nchar(synth_modules_name)))
    
    message('Lambda noise: ',lambda_noise)
    
    #read synths
    synth_module <- synth_modules[[synth_modules_name]]
    
    #read lognorm
    lognorm_est_counts <- synth_module$expmat
    
    #for linker
    regulator_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_module$regs)
    target_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_module$targs)
    
    
    #LINKER
    try({linker_output_synth <- TraRe::LINKER_run(lognorm_est_counts = lognorm_est_counts,
                                                  target_filtered_idx = target_filtered_idx, 
                                                  regulator_filtered_idx = regulator_filtered_idx,
                                                  NrCores = 32, Nr_bootstraps = 5, NrModules=10,
                                                  link_mode = methodlink, graph_mode = methodlink, FDR=FDR, onlymods=TRUE)
    
    
    saveRDS(linker_output_synth,paste0(getwd(),'/input_data/varying_noise_0.01_to_100_individuals/', methodlink,
                                       '_noise_variance_',lambda_noise,'.rds'))})
    
  }
  
}

file_samples <- paste0(getwd(),'/input_data/synth_modules_0.01_to_100.rds')
synth_modules <- readRDS(file_samples)
#call LINKER
difmodels_difnoise(synth_modules=synth_modules,methodlink_v = c('VBSR','LASSOmin','LM'),FDR=0.05)


difmodels_difsamples <- function(synth_modules, methodlink_v='LASSOmin',FDR=0.05){
  
  for (methodlink in methodlink_v){
    
    message('Working with.. ',methodlink)
    
    for (synth_modules_name in names(synth_modules)){
      
      nr_samples <- as.numeric(substr(synth_modules_name,9,nchar(synth_modules_name)))
      
      message('Nr samples: ',nr_samples)
      
      #read synths
      synth_module <- synth_modules[[synth_modules_name]]
      
      #read lognorm
      lognorm_est_counts <- synth_module$expmat
      
      #for linker
      regulator_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_module$regs)
      target_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_module$targs)
      
      
      #LINKER
      try({linker_output_synth <- TraRe::LINKER_run(lognorm_est_counts = lognorm_est_counts,
                                                    target_filtered_idx = target_filtered_idx, 
                                                    regulator_filtered_idx = regulator_filtered_idx,
                                                    NrCores = 64, Nr_bootstraps = 5, NrModules=10,
                                                    link_mode = methodlink, graph_mode = methodlink, FDR=FDR, 
                                                    onlymods=TRUE,only_train = TRUE)
      
      dir.create(paste0(getwd(),'/input_data/vary_samples/',methodlink))
      saveRDS(linker_output_synth,paste0(getwd(),'/input_data/vary_samples/', methodlink,'/',
                                         methodlink,'_nrsamples_',nr_samples,'.rds'))})
      
    }
  }
}


#Retrieve the samples file
file_samples <- paste0(getwd(),'/input_data/vary_samples/synth_modules_18_to_46.rds')
synth_modules <- readRDS(file_samples)
#call LINKER
difmodels_difsamples(synth_modules=synth_modules,methodlink_v = c('VBSR','LASSOmin','LM'),FDR=0.05)





