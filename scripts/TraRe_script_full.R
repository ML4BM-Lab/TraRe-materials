# SCRIPT TO RUN ALL STEPS OF TraRe on PROMOTE data

# OUTPUTS:
  # TraRe's GRN inference output for 10 and 50 bootstraps as 'linker_output.rds' files
  # Rewiring results for 10 bootstraps (directory with corresponding rewiring files)
  # Rewiring results for 50 bootstraps (only rewired sub-modules as 'fs_0.05.txt' file)
  # Rewiring at the gene level: Manuscript Table 1
  # Rewiring at the regulon level: '50b_regulons.rds' file and 'rewired_regulons_50b.xlsx' file for Supplementary material (extra)
   

# LIBRARIES ---------------------------------------------------------------

library(TraRe)
library(openxlsx)

# INPUT  -------------------------------------------------------------------
## Paths:
    wpath <- paste0(getwd(),'/TraRe-materials/')
    setwd(wpath)
    # PROMOTE
    # Expression matrix
    exp_p <- './input_data/mrm_sm_3500.txt'

    # Clinical info (Responder - non-responders)
    clinic_p <- './input_data/promote_v1_clinical.txt'

    # Gene info: drivers and targets
    gene_info_p <-'./input_data/promote_v1.gene_info.txt'


## Load data
    lognorm_est_counts <- read.table(exp_p)

    gene_info <- read.delim(gene_info_p)
    regulator_name <- gene_info[which(gene_info$regulator==1),1]
    target_name <- gene_info[which(gene_info$regulator==0),1]

    regulator_filtered_idx<-which(rownames(lognorm_est_counts)%in%regulator_name)
    target_filtered_idx <- which(rownames(lognorm_est_counts)%in%target_name)


# GRN inference method ----------------------------------------------------
  # 10 bootstraps
    linkeroutput <- LINKER_run(lognorm_est_counts =  as.matrix(lognorm_est_counts),
                                target_filtered_idx = target_filtered_idx,
                                regulator_filtered_idx = regulator_filtered_idx,
                                link_mode = "VBSR", #phase 1
                                graph_mode = "VBSR", #phase 2
                                NrModules = 100, #default
                                Nr_bootstraps = 10, #default
                                NrCores = 32)
    
    saveRDS(linkeroutput,file = './output/sm_3500_linker_output_10b.rds')
  
  # 50 bootstraps
    linkeroutput <- LINKER_run(lognorm_est_counts =  as.matrix(lognorm_est_counts),
                               target_filtered_idx = target_filtered_idx,
                               regulator_filtered_idx = regulator_filtered_idx,
                               link_mode = "VBSR", #phase 1
                               graph_mode = "VBSR", #phase 2
                               NrModules = 100, #default
                               Nr_bootstraps = 50, #default
                               NrCores = 32)
    
    saveRDS(linkeroutput,file = './output/sm_3500_linker_output_50b.rds')
    
# Rewiring ----------------------------------------------------------------
  # 10 Bootstraps
  prew_output <- preparerewiring(name = "rewiring_promote_sm3500_10b",
                               linker_output_p = './output/sm_3500_linker_output_10b.rds',
                               lognorm_est_counts_p = exp_p,
                               gene_info_p = gene_info_p,
                               phenotype_p = clinic_p,
                               outdir = './rewirings/', 
                               final_signif_thresh = 0.05,
                               nrcores = 4, 
                               use_graphs = T)
  runrewiring(prew_output)

  # At the gene level (50 bootstraps)
  impgenes <-  rewiring_gene_level(linker_output_p = './output/sm_3500_linker_output_50b.rds',
                                   lognorm_est_counts_p = exp_p,
                                   gene_info_p = gene_info_p,
                                   phenotype_p = clinic_p,
                                   fpath = "./output/50_bootstraps_rewmods.txt",
                                   final_sig_th = 0.05,
                                   include_cliques = FALSE,
                                   ImpTH = 0.05,
                                   nrcores = 4)
  # Save (Table 1)
  write.table(impgenes,file= './output/tables/table_rewired_tfs.tsv', row.names = TRUE, sep = "\t")
  
  # At the regulon level (50 bootstraps)
  sigmodules_p <- './output/50_bootstraps_rewmods_fs_0.05.txt'
  
  regulons <- rewiring_regulon_level(linker_output_p,
                                     lognorm_est_counts_p,
                                     gene_info_p,
                                     phenotype_p,
                                     sigmodules_p,
                                     rewired = TRUE,
                                     final_signif_thresh = 0.05)
  # Save output
  saveRDS(regulons,file = './output/sm_3500_50b_regulons.rds')
  
################################ EOF #############################################
                                     

