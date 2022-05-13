#Load linkeroutput

#pp <- 'G:/My Drive/PhD/PhDProjects(Papers)/TraRe/Trare - Main/Huanyao/'
# 
#linkeroutput <- readRDS(paste0(pp,'Module 3500/sm_3500_linker_output_50bootstraps.rds'))
# 
# #we need geneinfo as well
# 
#geneinfo <- read.delim(paste0(pp,'/Module 3500/promote_v1.gene_info.txt'))[,c('uniq_isos','regulator')]
# 
#rownames(geneinfo) <- geneinfo[,1]


graph_to_modules <- function(linkeroutput,geneinfo){
  
  #Ensure that rownames of geneinfo is the same as the uniq_isos
  rownames(geneinfo) <- geneinfo[,1]
  
  ## The structure we want to get is
  ## linkeroutput$modules[[link_mode]][[graph_mode]][[num_module]]$(target or regs)
  
  #For each module method:
  linkeroutput_new <- lapply(names(linkeroutput$modules),function(x){
    
    graph_modes <- names(linkeroutput$graphs[[x]])
    
    ## Select VBSR if more than 1 is available
    
    if (length(graph_modes)>1 & 'VBSR'%in%graph_modes){
      
        selected <- 'VBSR'
        
    }else{
      
      #Select the first one if there is only one or VBSR is not available
      selected <- graph_modes[1]
      
    }

    #Access to the selected graphs
    module_list <- lapply(seq_along(linkeroutput$graphs[[x]][[selected]]),function(y){
      
      #Save the graph
      graph <- linkeroutput$graphs[[x]][[selected]][[y]]
      
      #Select all the genes from the graph
      totgenes <- unique(names(igraph::V(graph)))
      
      #Select regulators
      regulators <- intersect(totgenes,rownames(geneinfo)[geneinfo[,2]==1])
      #Select target genes
      target_genes <- intersect(totgenes, rownames(geneinfo)[geneinfo[,2]==0])
      
      if (identical(regulators,character(0)) | identical(target_genes,character(0))){
        
        message('Module number ',y,' has been deleted')
        return(NULL)
        
      }
      
      list(regulators = regulators,
           target_genes = target_genes,
           bootstrap_idx = linkeroutput$modules[[x]][[y]]$bootstrap_idx)
      
    })
    
    #Create a copy containing the original indexes
    old2new <- sapply(module_list,function(x) !is.null(x))
    old2new <- seq_along(old2new)[old2new]
    
    #List result and name it
    module_list <- list(module_list[old2new])
    names(module_list) <- paste0(selected,'_modules')
    
    #Append old2new indexes
    module_list[[paste(selected,'original_index',sep='_')]] <- old2new

    #Drop empty graphs
    module_list[[selected]] <- module_list[[selected]][old2new]
    
    #Add new graphs
    module_list[[paste0(selected,'_graphs')]] <- linkeroutput[['graphs']][[x]][[selected]][old2new]
    
    return(module_list)
    
    
  })
    
    #Name them with the module methods
    names(linkeroutput_new) <- names(linkeroutput$modules)
    
    return(linkeroutput_new)
  
}

#ff <- graph_to_modules(linkeroutput,geneinfo)

