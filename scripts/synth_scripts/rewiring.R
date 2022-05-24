#Lets run rewiring
#load synth data
# synth_modules <- readRDS(paste0(getwd(),'/input_data/synth_modules_rewired_noise_1.2.rds'))
# #Check individual rewired modules
# #check_rewiring_linker_homemade(lognorm_est_counts,prep_sim$clinic,linker_output_synth,1)
# #Save lognorm to use in preparerewiring
# write.table(synth_modules$noise_1.2$expmat, 
#             paste0(getwd(),'/input_data/rewiring/synth_data_rew.txt'),
#             sep='\t',quote=F)

synth_data_p <- paste0(getwd(),'/input_data/rewiring/synth_data_rew.txt')
lognorm_est_counts <- as.matrix(read.delim(synth_data_p))
#First prepare

prep_noise <- '1.2'
for (prep_mode in c('LASSOmin','VBSR','LM')){
  
  prep_object <- TraRe::preparerewiring(name=paste0(c('synthdata',prep_mode),collapse='_'),
                                        linker_output_p = paste0(getwd(),'/input_data/models_noise_variance_1.2/',prep_mode,'_noise_variance_',prep_noise,'.rds'),
                                        lognorm_est_counts_p = paste0(getwd(),'/input_data/rewiring/synth_data_rew.txt'),
                                        gene_info_p = paste0(getwd(),'/input_data/rewiring/geneinfo.txt'),
                                        phenotype_p = paste0(getwd(),'/input_data/promote_clinical_ours_fixed.txt'),
                                        final_signif_thresh = 0.01,
                                        outdir = paste0(getwd(),'/input_data/rewiring/'),
                                        use_graphs = FALSE,
                                        nrcores=4)
  
  #check with the dave's method
  #check_rewiring_linker_dave(t(synth_modules$modules[[8]]),prep_sim$clinic)
  #Then run
  TraRe::runrewiring(prep_object)
}
