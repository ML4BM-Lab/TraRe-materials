#ROC as a function of pval
#Infer modules by using LINKER
generate_LINKER_pval <- function(synth_modules,method){
  
  lognorm_est_counts <- synth_modules$expmat
  #
  regulator_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_modules$regs)
  target_filtered_idx <- which(rownames(lognorm_est_counts)%in%synth_modules$targs)
  
  if (method=='LASSOmin'){
    
    for (Lambda in seq(2,8,1)){
      
      print(Lambda)
    
      #filepath
      filepath <- paste0(getwd(),'/input_data/vary_pvals/',method,'_lambda_',Lambda,'.rds')
      
      #Lets run LINKER
      linker_output_synth <- TraRe::LINKER_run(lognorm_est_counts = lognorm_est_counts,
                                               target_filtered_idx = target_filtered_idx, 
                                               regulator_filtered_idx = regulator_filtered_idx,
                                               NrCores = 32, Nr_bootstraps = 5, NrModules=10,
                                               link_mode = method, graph_mode = method, FDR = Lambda, onlymods=TRUE)
      
      saveRDS(linker_output_synth,filepath)
      
      
    }
  }else{ #VBSR,LM
    
    upperlimit <- length(synth_modules$regs)*length(synth_modules$targs)
    
    for (FDR in round(seq(upperlimit/(5*10),upperlimit/10,upperlimit/(5*10)))){
      
      print(FDR/upperlimit)
      
      #filepath
      filepath <- paste0(getwd(),'/try1/vary_pvals/',method,'_fdr_',FDR/upperlimit,'.rds')
      
      #Lets run LINKER
      linker_output_synth <- TraRe::LINKER_run(lognorm_est_counts = lognorm_est_counts,
                                               target_filtered_idx = target_filtered_idx, 
                                               regulator_filtered_idx = regulator_filtered_idx,
                                               NrCores = 32, Nr_bootstraps = 5, NrModules=10,
                                               link_mode = method, graph_mode = method, FDR = FDR, onlymods=TRUE)
      
      saveRDS(linker_output_synth,filepath)
      
    }
    
  }
  
}

#Evaluate LINKER
eval_linker <- function(linker_output, synth_mods, NrModules=10, TH=0.5, sankey=FALSE, quadgraphTH=0){
  
  #Retrieve lognorm 
  lognorm_counts <- synth_mods$expmat
  
  universe_size <- nrow(lognorm_counts)
  ref_drivers <- intersect(rownames(lognorm_counts),colnames(synth_mods$regprograms))
  
  #enrichment test
  hyptest <- function(mod1genes,mod2genes,universe_size){
    
    contig_tbl <- as.table(matrix(c(length(intersect(mod1genes, mod2genes)),
                                    length(setdiff(mod1genes, mod2genes)),
                                    length(setdiff(mod2genes, mod1genes)),
                                    universe_size - length(mod2genes) - length(mod1genes) + length(intersect(mod1genes, mod2genes))),
                                  ncol = 2, byrow = TRUE))
    #Show pvalue
    stats::fisher.test(contig_tbl, alternative = "g")$p.value
    
  }
  
  #jaccard index
  jaccindex <- function(mod1genes,mod2genes){
    
    return (length(intersect(mod1genes,mod2genes))/length(union(mod1genes,mod2genes)))
    
  }
  
  #calculate mode
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  #sankey diagram
  sank_diag <- function(adjmatrix,adjmat_link_gclustered){
    
    rownames_l <- length(rownames(adjmatrix))
    colnames_l <- length(colnames(adjmatrix))
    
    #Initialize Links
    Links <- data.frame()
    
    #For each Synth module
    for (i in seq(nrow(adjmatrix))){
      
      for (j in which(adjmatrix[i,]!=0)){
        Links <- rbind(Links,c(i-1,rownames_l+j-1,adjmatrix[i,j]*100))
      }
      
    }
    #Get infered module names
    #inf_names <- unique(colnames(adjmatrix)[Links[,2]-rownames_l+1])
    links_source <- Links[,1]+1
    
    #For each graph generated
    for (k in seq_along(adjmat_link_gclustered)){
      
      
      for (mod in adjmat_link_gclustered[[k]]){
        
        #there is no weight so we use 100
        Links <- rbind(Links,c(rownames_l+mod-1,rownames_l+colnames_l+k-1,100))
        
      }
      
    }
    
    colnames(Links) <- c('source','target','value')
    
    #Get nodes
    Nodes <- as.data.frame(c(rownames(adjmatrix),
                             paste0('I',seq_len(ncol(adjmatrix))),
                             paste0('G',seq_along(adjmat_link_gclustered))))
    colnames(Nodes) <- c('name')
    
    #Add group column
    Links <- cbind(Links,c(Nodes[links_source,],Nodes[Links[(length(links_source)+1):nrow(Links),2]+1,]))
    
    #Duplicate that column
    Links <- cbind(Links,Links[,4])
    
    #Rename
    colnames(Links) <- c('source','target','value','nodegroup','linkgroup')
    
    #modify group category in the {Gs} as a function of input Jaccard
    Gs_pos <- as.numeric(which(sapply(sapply(Links$nodegroup,function(x) grep('G',x)),length)==1))
    
    for (i in Gs_pos){
      #check which is the maximum Jaccard index entry 
      max_pos <- which(Links$target==Links[i,]$source)
      max_jacc <- Links[max_pos[which.max(Links[max_pos,'value'])],'nodegroup']
      Links[i,'linkgroup'] <- max_jacc
    }
    
    #now do the same with the nodes
    #duplicate column
    Nodes <- cbind(Nodes,Nodes)
    colnames(Nodes) <- c('name','colourname')
    
    Nodes_pos <- seq(nrow(adjmatrix)+ncol(adjmatrix)+1,nrow(Nodes))
    
    for (i in Nodes_pos){
      
      #check which source has the mode
      node_pos <- which(Links$nodegroup==Nodes[i,'colourname'])
      Nodes[i,'colourname'] <- getmode(Links$linkgroup[node_pos])
      
    }
    
    #define color
    #synth modules
    #iwanthue(5, 0, 360, 36, 180, 13, 73, plot=TRUE)  # intense
    synth_cols <- hues::iwanthue(nrow(adjmatrix),hmin=0, hmax=360, cmin = 36, cmax= 180, lmin=13, lmax=73)
    synth_linker_cols <- hues::iwanthue(ncol(adjmatrix),cmin = 20, cmax= 60, lmin=50, lmax=95)
    #transform this into a scale
    domain_string <- paste0(c("['",paste0(unique(Nodes$colourname),collapse="','"),"']"),collapse='')
    mycolor_domain <- paste0('.domain(',domain_string,')')
    
    range_string <- paste0(c("['",paste0(c(synth_cols,synth_linker_cols),collapse="','"),"']"),collapse='')
    mycolor_range <- paste0('.range(',range_string,')',collapse='') 
    
    mycolor <- paste0(c('d3.scaleOrdinal()',mycolor_domain,mycolor_range),collapse=' ')
    
    
    # Plot
    networkD3::sankeyNetwork(Links = Links, Nodes = Nodes, Source = "source",
                             Target = "target", Value = "value", NodeID = "name",
                             units = "Jaccard Index", fontSize = 12, nodeWidth = 20, iterations=150,
                             nodePadding = 5,LinkGroup = "linkgroup", NodeGroup='colourname',colourScale = mycolor,
                             sinksRight = FALSE, margin=0, fontFamily = "Arial")
    
  }
  
  #convert to graph
  #maybe its not a good idea to work with graphs as 
  #sizes are really small.
  #linkeroutput_g <- graph_to_modules(linker_output,geneinfo)
  
  #lets check if we have high enrichment on the modules
  #for each module in our linkeroutput
  
  #adjmatrix <- matrix(0,nrow=NrModules,ncol=length(linker_output$modules$VBSR))
  
  #generate the adjacency matrix
  #Columns
  #Among synthetic and inferred modules
  adjmatrix <- sapply(linker_output$modules[[1]],function(module){
    
    regs <- module$regulators
    targs <- module$target_genes
    real_genes <- c(regs,targs)
    
    #Rows
    sapply(synth_mods$modules,function(synthmod){
      
      
      synth_genes <- rownames(synthmod)
      jaccindex(real_genes,synth_genes)
      
      
    })
    
  })
  
  #Inferred (I) and real modules (R)
  rownames(adjmatrix) <- paste0('R',seq_along(synth_mods$modules))
  colnames(adjmatrix) <- paste0('I',seq_along(linker_output$modules[[1]]))
  
  ##Sankey Diagram
  #Adjmatrix among LINKER modules
  if (sankey){
    
    adjmatrix_linker <- sapply(linker_output$modules[[1]],function(linkmod1){
      
      regs1 <- linkmod1$regulators
      targs1 <- linkmod1$target_genes
      
      #Rows
      sapply(linker_output$modules[[1]],function(linkmod2){
        
        regs2 <- linkmod2$regulators
        targs2 <- linkmod2$target_genes
        
        jaccindex(c(regs1,targs1),c(regs2,targs2))
        
        
      })
      
    })
    adjmatrix_linker[adjmatrix_linker<TH] <- 0
    
    adjmat_link_g <- igraph::graph_from_adjacency_matrix(adjmatrix_linker,diag = FALSE,
                                                         weighted=TRUE,mode='undirected')
    adjmat_link_gclustered <- igraph::cluster_edge_betweenness(adjmat_link_g)
    
    #Generate the sankey diagram
    return(sank_diag(adjmatrix,adjmat_link_gclustered))
  }
  
  ##Quadgraphs
  
  message('Evaluating: ',quadgraphTH)
  #select the maximum
  if (inherits(quadgraphTH,'character')){
    
    #get the max and the position in the matrix
    max_pos <- apply(adjmatrix,2,which.max)
    
    for (column in seq(ncol(adjmatrix))){
      adjmatrix[-max_pos[column],column] <- 0 
    }
    
  }else{
    #Apply a threshold to it (0.2)
    adjmatrix[adjmatrix<quadgraphTH] <- 0
  }
  
  #Calculate the cuadrants graph
  #Iterate over the inferred (LINKER) modules
  quad_graphs <- lapply(seq_len(ncol(adjmatrix)),function(i){
    
    #Get the synth modules positions that have a non-zero Jaccard Index
    synth_pos_v <- which(adjmatrix[,i]!=0)
    
    #retrieve inferred regprog (LINKER) (=1 module)
    infered_regprog <- linker_output$modules[[1]][[i]]$regulatory_program
    infered_regprog_nonzero <- infered_regprog[infered_regprog!=0]
    #infered_regprog_nonzero_s <- sort(infered_regprog_nonzero,decreasing=TRUE)
    
    #For each synth module that has non-zero Jaccard indices
    #Create the matrix of points
    point_matrix <- matrix(0,2,1)
    
    #retrieve synth regprog (>=1 modules)
    for (synth_pos in synth_pos_v){
      
      #retrieve synth regprog (=1 modules as is in the loop)
      synth_regprog <- synth_mods$regprograms[synth_pos,]
      synth_regprog_nonzero <- synth_regprog[synth_regprog!=0]
      
      point_matrix <- cbind(point_matrix,sapply(ref_drivers,function(driver){
        
        driver_in_synth <- driver%in%names(synth_regprog_nonzero)
        driver_in_inferred <- driver%in%names(infered_regprog_nonzero)
        #conf_mat <- matrix(0,1,2)
        
        if (driver_in_synth & driver_in_inferred){
          #TP
          return (matrix(c(synth_regprog_nonzero[driver],infered_regprog_nonzero[driver]),1,2))
          
        }else if (driver_in_synth & !driver_in_inferred){
          #FN
          return (matrix(c(synth_regprog_nonzero[driver],0),1,2))
          
        }else if (!driver_in_synth & driver_in_inferred){
          #FP
          return (matrix(c(0,infered_regprog_nonzero[driver]),1,2))
          
        }else{
          #TN
          return (matrix(0,1,2))
        }
        
      }))
      
    }
    return(point_matrix[,-1])
    
  })
  
  #concat the matrices
  quad_graphs_c <- do.call(cbind,quad_graphs)
  
  #Retrieve TP,FP,FN
  confmat <- apply(quad_graphs_c,2,function(x){
    
    if (x[1]!=0 & x[2]!=0) return('TP')
    #not generated but inferred
    else if (x[1]!=0 & x[2]==0) return('FN')
    #not inferred but generated
    else if (x[1]==0 & x[2]!=0) return('FP')
    #not inferred nor generated
    else if (x[1]==0 & x[2]==0) return ('TN')
  })
  
  tp_fp_fn <- table(confmat)
  
  #fill the table
  for (fill in c('TP','TN','FP','FN')){
    if (!fill%in%names(tp_fp_fn)){
      tp_fp_fn[fill] <- 0
    }
  }
  
  return(list(confmat = confmat, confmat_table=tp_fp_fn,quadgraph=quad_graphs_c,adjmatrix = adjmatrix))
  
  
}

#Compute the tpr and fpr to build the roc curves as a function of p-values
roc_curves_pvals <- function(quadgraphTH,synth_modules){
  roc_curves <- function(folder_p,synth_modules,quadgraphTH,mode){
    
    roc <- matrix(0,3,1)
    
    for (file in list.files(folder_p)){
      
      message(folder_p,file)
      
      if (mode!='lasso'){
        pval <- as.numeric(stringr::str_extract(file,'[0-9].[0-9]+'))
      }else{
        pval <- as.numeric(stringr::str_extract(file,'[0-9]'))
      }
      
      #read the element
      linkfile <- readRDS(paste0(folder_p,file))
      
      #Call eval
      evalobject <- eval_linker(linkfile,synth_mods=synth_modules,quadgraphTH=quadgraphTH)
      
      #compute the roc curve
      tpf_sensit <- round(evalobject$confmat_table["TP"]/(evalobject$confmat_table["TP"]+evalobject$confmat_table["FN"]),2)
      fpr <- round(evalobject$confmat_table["FP"]/(evalobject$confmat_table["FP"]+evalobject$confmat_table["TN"]),2)
      
      roc <- cbind(roc,c(tpf_sensit,fpr,pval))
    }
    rownames(roc) <- c('TPR','FPR','Pval')
    
    #order
    if (mode=='lasso'){
      #First lambda is the maximum
      roc[3,1] <- 10
      rownames(roc) <- c('TPR','FPR','Lambda')
      roc <- roc[,order(roc[3,],decreasing = TRUE)]
      #normalize pval
      roc[3,] <- roc[3,]/max(roc[3,])
    }else{
      #vbsr/lm
      roc <- roc[,order(roc[3,])]
    }
  
    return(roc)
  }
  
  #Generate ROC Curves varying pval
  rocs_vbsr <- roc_curves(folder_p=paste0(getwd(),'/input_data/vary_pvals/vbsr/'),
                          synth_modules,quadgraphTH=quadgraphTH,mode='vbsr')
  rocs_lm <- roc_curves(folder_p=paste0(getwd(),'/input_data/vary_pvals/lm/'),
                        synth_modules,quadgraphTH=quadgraphTH,mode='lm')
  rocs_lasso <- roc_curves(folder_p=paste0(getwd(),'/input_data/vary_pvals/lasso/'),
                           synth_modules,quadgraphTH=quadgraphTH,mode='lasso')
  
  
  huescolors <- hues::iwanthue(3,hmax=120,cmin=30,cmax=80,lmin=35,lmax=80)
  
  
  
  #concats rows
  df_rows <- rbind(t(rocs_vbsr),t(rocs_lm),t(rocs_lasso))
  df_rows <- cbind(df_rows,c(rep('VBSR',ncol(rocs_vbsr)),
                             rep('LM',ncol(rocs_lm)),
                             rep('LASSO',ncol(rocs_lasso))))
  df_rows <- cbind(df_rows,df_rows[,3])
  colnames(df_rows) <- c('TPR','FPR','Pval','Model','Lambda')
  
  #remove duplicities
  df_rows[which(!df_rows[,'Model']%in%c('VBSR','LM')),'Pval'] <- NA
  df_rows[which(df_rows[,'Model']%in%c('VBSR','LM')),'Lambda'] <- NA
  
  #filter by FPR
  FPR_th <- 0.1
  #df_rows_filtered <- df_rows[apply(df_rows,1,function(x) if (as.numeric(x[2])<=FPR_th){TRUE}else{FALSE}),]
  df_rows_filtered <- df_rows
  
  #add pval cond to remove 0
  df_rows_filtered <- cbind(df_rows_filtered,df_rows_filtered[,'Pval'])
  colnames(df_rows_filtered) <- c('TPR','FPR','Pval','Model','Lambda','PvalCond')
  
  pval <- as.numeric(df_rows_filtered[,'PvalCond'])
  df_rows_filtered[which(pval==0),'PvalCond'] <- NA
  
  #Define the range for the scale
  pvalsize <- as.numeric(na.omit(unique(as.numeric(df_rows_filtered[,'PvalCond']))))
  pvalsize_l <- c(min(pvalsize),mean(pvalsize),max(pvalsize))
  
  #add the lasso points
  lasso_p <- df_rows_filtered[,'Model']=='LASSO'
  lasso_df <- df_rows_filtered[lasso_p,]
  
  library(ggplot2)
  roc_curve_gg <- ggplot(as.data.frame(df_rows_filtered),aes(x=as.numeric(FPR), y=as.numeric(TPR),fill=Model)) + 
    
    #VBSR
    ggrepel::geom_label_repel(aes(label = Lambda),
                              size=4,
                              box.padding  =  1.2,
                              point.padding = 0.5,
                              segment.size = 0.5,
                              segment.color = 'grey50',
                              max.overlaps = 50,
                              na.rm = TRUE) +
    geom_point(aes(size=as.numeric(PvalCond)),na.rm = TRUE,color='black',shape=21,stroke=1)+
    geom_point(data=as.data.frame(lasso_df),na.rm = TRUE)+
    scale_size("p-values", range=c(8, 1), breaks=pvalsize)+
    #scale_size_continuous(range = c(1,5),breaks= -pvalsize_l, name = "Pvalue")+
    #geom_point() + 
    geom_line(aes(color=Model))+
    labs(x='FPR', y = 'TPR') +
    theme_classic() +
    scale_fill_manual(values = hues::iwanthue(length(unique(df_rows_filtered[,"Model"])),
                                              hmax=120,cmin = 30,cmax = 80,lmin=35,lmax=80),guide=FALSE) +
    scale_color_manual(values = hues::iwanthue(length(unique(df_rows_filtered[,"Model"])),
                                               hmax=120,cmin = 30,cmax = 80,lmin=35,lmax=80)) +
    ggtitle('ROC curves as a function of pvalues threshold')+
    theme(plot.title = element_text(hjust = 0.5))+
    xlim(0,0.2)
    #scale_x_continuous(breaks = seq(0,0.2,0.025))
  
  pdf(file=paste0(getwd(),'/output/figures/roc_pvals_',quadgraphTH,'.pdf'))
  print(roc_curve_gg)
  dev.off()
  
  
}

synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_0.01_to_100.rds'))

#Generate Pvals for noise lambda=1.2
#Define the method
# for (method in c('VBSR','LASSOmin','LM')){
#   print(method)
#   generate_LINKER_pval(synth_modules = synth_modules$noise_1.2, method=method)
#   
# }
#   

for (quadgraphTH in c(0,0.5,0.8,'max')){
  ## Generate ROC curves as a function of pvals
  df <- roc_curves_pvals(quadgraphTH=quadgraphTH,synth_modules=synth_modules$noise_1.2)
}


