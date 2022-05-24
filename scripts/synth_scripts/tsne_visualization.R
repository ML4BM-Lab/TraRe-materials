library(ggplot2)

plot_lognorm <- function(synth_modules,mode='t-SNE',run=FALSE,linker_output_synth=FALSE){
  
  #Lets plot using numbers as synth clusters and colors as LINKER inferred clusters
  labels <- list(colors=c(),numbers=c())
  
  if (run){
    
    for (i in c(2,1)){#1 for colors, 2 for numbers
      
      if (i==1){
        
        #Keep only genes from the specific run
        mods <- which(sapply(linker_output_synth$modules[[1]],function(x) x$bootstrap_idx)==run)
        palett <- hues::iwanthue(length(mods))
        
        #Targets/module
        gene_labels <- sapply(mods, function(i){
          linker_output_synth$modules[[1]][[i]]$target_genes
        })
        
      }else{
        
        mods <- seq_along(synth_modules$modules)
        gene_labels <- lapply(synth_modules$modules,rownames)
        
        
      }
      
      #Build lognorm_counts from the synth_modules$expmat
      lognorm_est_counts <- synth_modules$expmat[unique(unlist(gene_labels)),]
      
      #lognorm_target_rownames
      rnames_target_lognorm <- sort(setdiff(rownames(lognorm_est_counts),unique(unlist(synth_modules$drivers_mod))))
      
      #dict gene-cluster it belongs to
      gene_cluster_dict <- sapply(rnames_target_lognorm,
                                  function(x) which(sapply(seq(mods),
                                                           function(y) x%in%gene_labels[[y]])))
      
      #Compute colors and numbers
      if (i==1){
        labels[[i]] <- sapply(rnames_target_lognorm,function(x) palett[gene_cluster_dict[x]])
      }else{
        labels[[i]] <- sapply(rnames_target_lognorm,function(x) gene_cluster_dict[x])
      }
    }
    
  }else{
    
    #We will be working only with target genes as drivers can be repeated
    #Retrieve lognorm
    lognorm_est_counts <- synth_modules$expmat
    
    gene_labels <- lapply(synth_modules$modules,rownames)
    
    mods <- length(synth_modules$modules)
    palett <- hues::iwanthue(mods)

    #Build lognorm_counts from the synth_modules$expmat
    lognorm_est_counts <- synth_modules$expmat[unique(unlist(gene_labels)),]

    #lognorm_target_rownames
    rnames_target_lognorm <- sort(setdiff(rownames(lognorm_est_counts),unique(unlist(synth_modules$drivers_mod))))

    #dict gene-cluster it belongs to
    gene_cluster_dict <- unlist(sapply(rnames_target_lognorm,
                                function(x) which(sapply(seq(mods),
                                                         function(y) x%in%gene_labels[[y]]))))
    
    labels <- sapply(rnames_target_lognorm,function(x) palett[gene_cluster_dict[x]])
  }
  
  
  if (mode=='t-SNE'){
    
    if (run){
    #compute t-sne
    set.seed(42)
    try_tsne <- Rtsne::Rtsne(lognorm_est_counts[rnames_target_lognorm,])
    #try_tsne <- M3C::tsne(t(lognorm_est_counts))
    #try_tsne <- tsne::tsne(lognorm_est_counts,max_iter = 200)
    dataset <- cbind(as.data.frame(try_tsne$Y),labels$numbers)
    colnames(dataset) <- c('dim1','dim2','numbers')
    
    #Plot it
    return(ggplot2::ggplot(dataset, aes(x=dim1, y=dim2,label=numbers)) +
             geom_text(size=4,col=labels$colors)+
             labs(x= paste(mode,'1'), y= paste(mode,'2')) +
             theme_classic())
    
    }else{
      set.seed(42)
      try_tsne <- Rtsne::Rtsne(lognorm_est_counts[rnames_target_lognorm,])
      #try_tsne <- M3C::tsne(t(lognorm_est_counts))
      #try_tsne <- tsne::tsne(lognorm_est_counts,max_iter = 200)
      dataset <- as.data.frame(try_tsne$Y)
      colnames(dataset) <- c('dim1','dim2')
      
      #Plot it
      return(ggplot2::ggplot(dataset, aes(x=dim1, y=dim2)) +
               geom_point(size=2,shape=20,col=labels)+
               labs(x= paste(mode,'1'), y= paste(mode,'2')) +
               theme_classic())
    }
    
    
  }else if (mode=='PCA'){
    #Check with pca
    set.seed(42)
    decomp <- prcomp(lognorm_est_counts[rnames_target_lognorm,], scale = TRUE)
    
    dataset <- as.data.frame(decomp$x[,1:2])
    colnames(dataset) <- c('dim1','dim2')
    
    #Plot it
    return(ggplot2::ggplot(dataset, aes(x=dim1, y=dim2)) +
             geom_point(size=2,shape=20,col=labels)+
             labs(x= paste(mode,'1'), y= paste(mode,'2')) +
             theme_classic())
    
  }
  
  
}


#Plot synth_modules
synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_0.01_to_100.rds'))

#noise_0.4 & noise_1.2
for (noise in c('noise_0.4','noise_1.2')){
  
#Plot the PCA
pdf(paste0(getwd(),'/output/figures/synth_modules_pca_',noise,'.pdf'))
print(plot_lognorm(synth_modules[[noise]],mode='PCA'))
dev.off()

}
