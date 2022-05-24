#Evaluate LINKER + Sankey Diagram
eval_linker <- function(linker_output, synth_mods, NrModules=10, TH=0.5, sankey=FALSE, quadgraphTH=0, boxplot=FALSE){
  
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
  
  if (boxplot){
    return(adjmatrix)
  }
  
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

synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_rewired_noise_1.2.rds'))

for (method in c('LM','LASSOmin','VBSR')){

  linker_output_synth <- readRDS(paste0(getwd(),'/input_data/models_noise_variance_1.2/',method,'_noise_variance_1.2.rds'))
  #these have to be saved manually
  sankey_method <- eval_linker(linker_output_synth,synth_modules$noise_1.2,TH=0.5,sankey = TRUE)

}

#BOXPLOT
#Generate boxplot to see differences between synthetic and inferred modules among VBSR,LM and LASSO
boxplot_points_all_the_same_first_fifty <- function(){
  
  models <- c('LM','LASSOmin','VBSR')
  
  bp_v <- c()
  for (mode in models){
    bp_tmp <- c()
    linker_output_synth <- readRDS(paste0(getwd(),'/input_data/models_noise_variance_1.2/',mode,'_noise_variance_1.2.rds'))
    for (i in seq(10)){
      bp_tmp <- c(bp_tmp,eval_linker(linker_output_synth,synth_modules$noise_1.2,TH=0.5,sankey = TRUE,boxplot=TRUE)[i,])
    }
    bp_v <- c(bp_v,bp_tmp)
  }
  
  mode_v <- c(rep(models[1],36*10),rep(models[2],50*10),rep(models[3],50*10))
  rep_v <- function(j){
    return(unlist(lapply(seq(10),function(i) rep(i,j))))
  }
  
  #remove 0s
  #bp_v[bp_v==0] <- NA
  bp_mat <- cbind(bp_v,mode_v)
  get_fifties <- function(mat){
    
    new_mat <- matrix(NA,1,2)
    
    for (mode in unique(mat[,2])){
      
      pos <- mat[,2]==mode
      v <- as.numeric(mat[,1][pos])
      
      new_v <- sort(v,decreasing = TRUE)[1:50]
      new_v_mode <- cbind(new_v,rep(mode,length(new_v)))
      
      new_mat <- rbind(new_mat,new_v_mode)
      
      
    }
    
    colnames(new_mat) <- c('Boxplot','Mode')
    return(new_mat[-1,])
    
  }
  
  matt <- get_fifties(bp_mat)
  
  plot_box <- function(matt){
    
    colors <- hues::iwanthue(n=3,hmin = 181,hmax=285,cmin=40,cmax=75,lmin=48,lmax=80)
    colorsv2 <- hues::iwanthue(n=3,hmin = 143,hmax=238,cmin=40,cmax=78,lmin=45,lmax=78)
    
    my_comparisons <- list( c("LM", "LASSOmin"), c("LM", "VBSR"),c("LASSOmin", "VBSR"))
    ggplot(as.data.frame(matt), aes(x= factor(Mode,levels=models),
                                    y=as.numeric(Boxplot))) +
      
      geom_boxplot(fill=colors) +
      ggsignif::geom_signif(comparisons = list( c("LM", "LASSOmin")),test='wilcox.test',
                            map_signif_level = TRUE, textsize = 6,
                            vjust = .5, y_position = 1.05) +
      
      ggsignif::geom_signif(comparisons = list( c("LM", "VBSR")),test='wilcox.test',
                            map_signif_level = TRUE, textsize = 6,
                            y_position = 1.15, vjust= .5) +
      
      ggsignif::geom_signif(comparisons = list( c("LASSOmin", "VBSR")),test='wilcox.test',
                            map_signif_level = TRUE, textsize = 6, vjust= .5) +
      #ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") + 
      labs(x='Models', y = 'Jaccard Index') +
      theme_classic() 
  }
  
  pdf(paste0(getwd(),'/output/figures/sankey_jaccard_index_models.pdf'))
  print(plot_box(matt))
  dev.off()
  
  
}
#boxplot + significance
boxplot_points_all_the_same_first_fifty()
