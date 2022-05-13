# SCRIPT TO GENERATE SUB-MODULES WITH DIFFERENT PARAMETER IN NUMBER OF MODULES PER BOOTSTRAP.
# AND THE CORESPONDING PLOT FOR THE SUPPLEMENTARY FIGURE 1

# LIBRARIES ---------------------------------------------------------------
library(TraRe)
library(dplyr)
library(ggplot2)
library(grpubr)
# INPUT  -------------------------------------------------------------------
  ## Paths:
  wpath <- paste0(getwd(),'/TraRe_NAR/')
  setwd(wpath)
  # PROMOTE
  # Expression matrix
  exp_p <- './input_data/mrm_sm_3500.txt'
  
  # Clinical info (Responder - non-responders)
  clinic_p <- '/input_data/promote_clinical_ours_improved.txt'
  
  # Gene info: drivers and targets
  gene_info_p <-'./input_data/promote_v1.gene_info.txt'
  
  
  ## Load data
  lognorm_est_counts <- read.table(exp_p)
  
  gene_info <- read.delim(gene_info_p)
  regulator_name <- gene_info[which(gene_info$regulator==1),1]
  target_name <- gene_info[which(gene_info$regulator==0),1]
  
  regulator_filtered_idx<-which(rownames(lognorm_est_counts)%in%regulator_name)
  target_filtered_idx <- which(rownames(lognorm_est_counts)%in%target_name)
  
  
# GRN inference loop ----------------------------------------------------
  nmods <- c(50,100,200,300,500)
  for (n in nmods){
  linkeroutput <- LINKER_run(lognorm_est_counts =  lognorm_est_counts,
                             target_filtered_idx = target_filtered_idx,
                             regulator_filtered_idx = regulator_filtered_idx,
                             link_mode = "VBSR", #phase 1
                             graph_mode = "VBSR", #phase 2
                             NrModules = n, #default
                             Nr_bootstraps = 5, #default
                             corrClustNrIter = 500,
                             NrCores = 32)
  
  saveRDS(linkeroutput, './output/linker_output_5b_',n,'mods','.rds')
  }


# PLOT --------------------------------------------------------------------

  #### LOAD DATA AND EXTRACT STATS ####
  
  lo50 <- readRDS('.output/linker_output_5b_50mods.rds')
  lo200 <- readRDS('.output/linker_output_5b_100mods.rds')
  lo100 <- readRDS('.output/linker_output_5b_200mods.rds')
  lo300 <- readRDS('.output/linker_output_5b_300mods.rds')
  lo500 <- readRDS('.output/linker_output_5b_500mods.rds')
  
  lo_list <- list(lo50,lo100,lo200, lo300, lo500)
  
  get_stats <- function(lo_list){
    out <- lapply(seq_along(lo_list), function(lo){
      x <- lo_list[[lo]]
      NrModules <- sapply(x$raw_results$VBSR$bootstrapResults, function(y) y$NrModules)
      genes <- lapply(x$modules$VBSR, function(m){
        n_regs_m <- length(m$regulators)
        n_targs_m <- length(m$target_genes)
        return(list(n_regs_m=n_regs_m,n_targs_m=n_targs_m))
      })
      genes_df <- as.data.frame(matrix(unlist(genes),ncol=2, byrow=T))
      colnames(genes_df) <- c("regulators", "targets")
      genes_df$bootstrap <- as.factor(unlist(sapply(seq(5),function(x) rep(x,NrModules[x]))))
      genes_df$group <- rep(lo,dim(genes_df)[1])
      return(genes_df)
    })  
    return(do.call(rbind,out)) 
  }
  
  stats_mods <- get_stats(lo_list)
  stats_mods$group <-factor(stats_mods$group, labels=c("50","100", "200", "300", "500" ))
  
  ##### PLOT #####
  p1 <- ggplot(stats_mods, aes(x=group, fill=group))+
    geom_bar(stat="count",show.legend = FALSE)+
    theme_bw(base_size = 25)+
    theme(plot.title = element_text(hjust = 0.5,size = 25, family="sans"),
          aspect.ratio = 1,
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    ylab('')+
    xlab('')+
    ggtitle('Number of modules')

  # Regulons
  p2 <- ggplot(stats_mods, aes(x = group))+
    theme_bw(base_size = 25)+
    theme(plot.title = element_text(hjust = 0.5,size = 25, family="sans"),
          aspect.ratio = 1,
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    geom_jitter(aes(y=regulators,colour=group, fill=group),show.legend = FALSE)+
    ylab('')+
    xlab('')+
    ggtitle('Transcription Factors')
  
  
  
  # Targets
  p3 <-ggplot(stats_mods, aes(x= group ))+
    theme_bw(base_size = 25)+
    theme(plot.title = element_text(hjust = 0.5,size = 25, family="sans"),
          legend.title = element_text(size=13),
          legend.text = element_text(size=11),aspect.ratio = 1,
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    geom_violin(aes(y=targets,colour=group, fill=group))+
    ylab('')+
    xlab('')+
    ggtitle('Target genes')+
    scale_fill_discrete(name = "Number of\npre-set modules",labels=c("50 mods", "100 mods", "200 mods", "300 mods", "500 mods"))+
    scale_colour_discrete(name = "Number of\npre-set modules", labels=c("50 mods", "100 mods", "200 mods", "300 mods", "500 mods"))
   
  
  pdf('./output/suppl_fig_griding_over_nrmods.pdf', width=16,height = 6,onefile=FALSE) 
  ggpubr::ggarrange(p1,p2,p3,nrow=1,ncol=3, common.legend = T, legend = "right",align = "hv")+
    annotate('text', x = 0.05, y = 0.9,label = 'A', size = 12)+
    annotate('text', x = 0.35, y = 0.9,label = 'B', size = 12)+
    annotate('text', x = 0.65, y = 0.9,label = 'C', size = 12)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  dev.off()
################################ EOF #############################################