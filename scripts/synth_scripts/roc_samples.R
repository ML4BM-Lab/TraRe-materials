#ROC as a function of samples

roc_curves_samples <- function(synth_modules,folder_p,quadgraphTH_v='max'){
  
  for (quadgraphTH in quadgraphTH_v){
    
    if (quadgraphTH!='max'){
      quadgraphTH <- as.numeric(quadgraphTH)
    }
    #message('Evaluating threshold: ',quadgraphTH)
    
    compute_tpr_fpr <- function(eval_object){
      tpf_sensit <- round(eval_object$confmat_table["TP"]/(eval_object$confmat_table["TP"]+eval_object$confmat_table["FN"]),5)
      fpr <- round(eval_object$confmat_table["FP"]/(eval_object$confmat_table["FP"]+eval_object$confmat_table["TN"]),5)
      return(c(tpf_sensit,fpr))
    }
    
    #Initialize ROC_mat
    ROC_mat <- matrix(NA,length(names(synth_modules)),6) #TPR,FPR for VBSR,LASSO and LM
    
    #assign colnames y rownames
    colnames(ROC_mat) <- c('VBSR_TPR','VBSR_FPR',
                           'LM_TPR','LM_FPR',
                           'LASSOmin_TPR','LASSOmin_FPR')
    rownames(ROC_mat) <- names(synth_modules)
    
    for (nrsamples in names(synth_modules)){
      
      message('samples ',nrsamples)
      #message('Evaluating ',nrsamples)
      synth_module <- synth_modules[[nrsamples]]
      
      for (method in c('VBSR','LM','LASSOmin')){
        
        message('Evaluating ',method)
        
        #Check if file exists
        file_p <- paste0(folder_p,method,'/',method,'_nr',nrsamples,'.rds')
        
        if (file.exists(file_p)){
          
          linker_output_synth <- readRDS(file_p)
          lognorm_est_counts <- synth_module$expmat
          
          eval_object <- eval_linker(linker_output = linker_output_synth,synth_mods = synth_module, 
                                     TH=0.50,sankey=FALSE,quadgraphTH=quadgraphTH)
          
          roc_params <- compute_tpr_fpr(eval_object)
          #message(roc_params)
          
          #Fill the matrix
          ROC_mat[nrsamples,paste(method,'TPR',sep='_')] <- roc_params[1]
          ROC_mat[nrsamples,paste(method,'FPR',sep='_')] <- roc_params[2]
        }
        
      }
      
    }
    
    #Format
    ROC_format <- matrix(NA,nrow(ROC_mat)*3,4) #TPR,FPR for VBSR,LASSO and LM
    colnames(ROC_format) <- c('TPR','FPR','Method','Samples')
    #Fill it
    for (method in c('VBSR','LM','LASSOmin')){
      i<-which(c('VBSR','LM','LASSOmin')%in%method)
      
      ROC_format[seq(nrow(ROC_mat)*(i-1)+1,nrow(ROC_mat)*i),'TPR'] <- ROC_mat[,paste0(c(method,'TPR'),collapse='_')]
      ROC_format[seq(nrow(ROC_mat)*(i-1)+1,nrow(ROC_mat)*i),'FPR'] <- ROC_mat[,paste0(c(method,'FPR'),collapse='_')]
      ROC_format[seq(nrow(ROC_mat)*(i-1)+1,nrow(ROC_mat)*i),'Method'] <- rep(method,nrow(ROC_mat))
      ROC_format[seq(nrow(ROC_mat)*(i-1)+1,nrow(ROC_mat)*i),'Samples'] <- substr(rownames(ROC_mat),9,nchar(rownames(ROC_mat)))
    }
    
    ##Plot
    library(ggplot2)
    scatter <- ggplot(as.data.frame(ROC_format), aes(x=as.numeric(FPR),
                                                     y=as.numeric(TPR),
                                                     fill=Method,
                                                     label=as.numeric(Samples))) +
      
      geom_point(na.rm = TRUE,stroke=1.2,shape=21,size=5) +
      geom_text(na.rm=TRUE,aes(label=as.numeric(Samples)),hjust=-0.5, vjust=-0.5)+
      
      # scale_size("Nº-samples", range=c(0.1, 10),
      #           breaks=as.numeric(unique(ROC_format[,'Samples'])),
      #           labels=unique(ROC_format[,'Samples']))+
      
      labs(x='FPR', y = 'TPR') +
      theme_classic() +
      
      ggtitle(paste0('ROC curves as a function of samples - Noise: 1.2 -  ','Threshold: ',quadgraphTH))+
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(name='Method',values = hues::iwanthue(3, hmax=120,cmin = 30,cmax = 80,lmin=35,lmax=80))+
      guides(fill = guide_legend(override.aes = list(shape = 21)), size = guide_legend(override.aes = list(shape = 21)))
    
    pdf(file=paste0(getwd(),'/output/figures/roc_samples_',quadgraphTH,'.pdf'))
    print(scatter)
    dev.off()
  }
  
}

#Generate plots
synth_modules <- readRDS(paste0(getwd(),'/input_data/vary_samples/synth_modules_18_to_46_samples.rds'))
#Generate ROC curves
folder_p<- paste0(getwd(),'/input_data/vary_samples/')

#evaluate multiple thresholds
for (quadgraphTH in c(0,0.5,0.8,'max')){
roc_curves_samples(synth_modules=synth_modules,folder_p=folder_p,quadgraphTH_v=quadgraphTH)
}
