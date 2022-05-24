#ROC as a function of noise
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
    
    
    saveRDS(linker_output_synth,paste0(getwd(),'/input_data/varying_noise_0.01_to_100_individual/', methodlink,
                                       '_noise_variance_',lambda_noise,'.rds'))})
    
  }
  
}
roc_curves_noise <- function(folder_p,work_seq,synth_modules,quadgraphTH_v){
  
  dfs_lapply <- lapply(quadgraphTH_v,function(quadgTH){
    
    if (!quadgTH=='max'){
      quadgTH <- as.numeric(quadgTH)
    }
    
    message('Evaluating ',quadgTH,' noise')
    
    #Initialize
    folder_path <- paste0(getwd(),folder_p)
    working_seq <- work_seq
    
    ROC_mat <- matrix(NA,length(working_seq),6) #TPR,FPR for VBSR,LASSO and LM
    
    #assign colnames y rownames
    colnames(ROC_mat) <- c('VBSR_TPR','VBSR_FPR',
                           'LM_TPR','LM_FPR',
                           'LASSO_TPR','LASSO_FPR')
    
    rownames(ROC_mat) <- as.character(working_seq)
    
    folder_path_files <- list.files(folder_path)
    
    for (lambda_noise in working_seq){
      
      message('Evaluating noise level ',lambda_noise)
      
      #select our synth module
      synth_module <- synth_modules[[paste0('noise_',lambda_noise)]]
      
      if (length(synth_module$modules)!=10){
        next
      }
      
      lognorm_est_counts <- synth_module$expmat
      
      for (model_i in seq(3)){
        
        model <- c('VBSR','LM','LASSOmin')[model_i]
        message('Model ', model)
        
        #grep the model and noise
        filep <- paste0(model,'_noise_variance_',lambda_noise,'.rds')
        
        #if the file exists
        if (file.exists(paste0(folder_path,filep))){
          
          #load linkeroutput
          linker_output_synth <- readRDS(paste0(folder_path,filep))
          
          #call evallinker
          eval_object <- eval_linker(linker_output=linker_output_synth,synth_mods=synth_module,
                                     TH=0.50,sankey=FALSE,quadgraphTH=quadgTH)
          
          #compute the roc curve
          tpf_sensit <- round(eval_object$confmat_table["TP"]/(eval_object$confmat_table["TP"]+eval_object$confmat_table["FN"]),2)
          fpr <- round(eval_object$confmat_table["FP"]/(eval_object$confmat_table["FP"]+eval_object$confmat_table["TN"]),2)
          
          row_i <- which(working_seq%in%lambda_noise)
          ROC_mat[row_i,c(model_i*2-1,model_i*2)] <- c(tpf_sensit,fpr)
          
        }
        
      }
      
      
    }
    
    
    df_lapply <- rbind(rbind(ROC_mat[,c(1,2)],ROC_mat[,c(3,4)]),ROC_mat[,c(5,6)])
    colnames(df_lapply) <- c('TPR','FPR')
    df_lapply <- cbind(df_lapply,c(rep('VBSR',length(work_seq)),
                                   rep('LM',length(work_seq)),
                                   rep('LASSOmin',length(work_seq))))
    #add noise
    df_lapply <- cbind(df_lapply,working_seq)
    
    #add method
    df_lapply <- cbind(df_lapply,rep(paste0('TH_',quadgTH),nrow(df_lapply)))
    colnames(df_lapply) <- c('TPR','FPR','method','noise','th_mode')
    
    return(df_lapply)
    
  })
  
  #Concat dfs
  concat_dfs <- do.call(rbind,dfs_lapply)
  
  #rename to df
  df <- concat_dfs
  
  #duplicate noise
  df <- cbind(df,df[,4])
  colnames(df) <- c('TPR','FPR','method','noise_low','th_mode','noise_large')
  
  #define noise breaks
  noise_th <- sapply(as.numeric(df[,'noise_low']),function(x) if (x<10){TRUE}else{FALSE})
  df[!noise_th,'noise_low'] <- as.character(max(as.numeric(df[noise_th,'noise_low'])))
  df[noise_th,'noise_large'] <- NA
  
  #noise_unique <- unique(as.numeric(df[,'noise']))
  noise_low_breaks <- as.numeric(quantile(na.omit(unique(as.numeric(df[,'noise_low'])))))
  noise_low_breaks_norm <- noise_low_breaks/max(noise_low_breaks)
  
  library(ggplot2)

  for (th in unique(df[,'th_mode'])){
    
    th_df <- df[df[,'th_mode']==th,]
    
    scatter <- ggplot(as.data.frame(th_df), aes(x=as.numeric(FPR),
                                                y=as.numeric(TPR),
                                                fill=method,
                                                size=as.numeric(noise_low)/max(as.numeric(noise_low)))) +
      
      geom_point(na.rm = TRUE,stroke=1.2,shape=21) +
      
      scale_size("\u03c3-noise", range=c(0.1, 6),breaks=noise_low_breaks_norm[-1],labels=noise_low_breaks[-1])+
      
      labs(x='FPR', y = 'TPR') +
      theme_classic() +
      
      ggrepel::geom_label_repel(aes(label = noise_large),
                                size=3,
                                box.padding   = 3,
                                point.padding = 0.5,
                                segment.color = 'grey50',
                                max.overlaps=500,
                                na.rm = TRUE,
                                direction='y',
                                show.legend = FALSE) +
      ggtitle('ROC curves as a function of noise variance within synthetic modules')+
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(name='Method',values = hues::iwanthue(length(unique(th_df[,3])),
                                                              hmax=120,cmin = 30,cmax = 80,lmin=35,lmax=80))+
      guides(fill = guide_legend(override.aes = list(shape = 21)), size = guide_legend(override.aes = list(shape = 21)))
    
      pdf(file=paste0(getwd(),'/output/figures/roc_noise_',th,'.pdf'))
      print(scatter)
      dev.off()
    
  }
  
  return(df)
  
  
}


synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_0.01_to_100.rds'))

#call for lasso,lm and vbr
# difmodels_difnoise(synth_modules, methodlink='LASSOmin')
# difmodels_difnoise(synth_modules, methodlink='VBSR')
# difmodels_difnoise(synth_modules, metodlink='LM')


## Generate ROC curves as a function of noise

#Generate ROC fixing pval (0.05 and lambdamin), varying noise variance
lambda_noise <- sort(as.numeric(sapply(names(synth_modules),function(x) substr(x,7,nchar(x)))))
ROC_noise_mat <- roc_curves_noise(folder_p='/input_data/varying_noise_0.01_to_100_individual/',
                                  work_seq=lambda_noise,synth_modules,quadgraphTH_v = c(0,0.5,0.8,'max'))

