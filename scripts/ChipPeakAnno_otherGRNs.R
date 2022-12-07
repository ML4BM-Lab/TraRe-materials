
### ------------------------------------------------------------
# library(ggplot2)
# library(ChIPpeakAnno)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(IRanges)
# library(xlsx)
# library(fmsb)

### ------LOAD CHIPPEAK DATA---------------------------------
## See CHIPSEQ.R script to see how these objects were generated
print("Loading data")
gencode33all <- readRDS('./DF_gencode.v33.annotation_filtered.rds')
G19AnnotChIP_RL <- readRDS('./GRanges.gencode.v33.annotation_filtered.rds')
all_peaks <- readRDS("./remap2022_chippeaks.rds")
all_TFs<- unique(all_peaks$peakName)

### FUNCTION FOR STATS-----------------------------------------------

get_stats_chipseq <- function(tf,ChipPeak_results,regulon_results,alternative="greater",out_targets=F){
        # Check
        ifelse(tf!=unique(ChipPeak_results[ChipPeak_results$peakName==tf,c("peakName")]), stop("Some misstakes in Chipseq data"),T)
        
        ifelse(tf%in%unique(regulon_results$TF),T,stop("Some misstakes in Chipseq data"))
        
        targets_in_ChIP <- unique(ChipPeak_results[ChipPeak_results$peakName==tf,c("gene_name")]) # SUCCESS
        all_regulon_targets <- regulon_results[regulon_results$TF==tf,]

        # Now select those where, having edges in both R and NR GRNs have the top 10% difference (it was already ordered)
        all_regulon_targets <- all_regulon_targets[order(all_regulon_targets$diff_importance, decreasing = T),]
        regulon_targets <- all_regulon_targets[1:round(nrow(all_regulon_targets)*0.25),2]
        
        # # PHYPER
        # N_hyper <- length(total_universe_targets) # UNIVERSE

        # n_hyper <- length(targets_in_ChIP) # SUCCESS

        # m_hyper <- N_hyper-n_hyper # NOT SUCCESS
 
        # k_hyper <- length(regulon_targets)# DRAWNS

        # q_hyper <- length(intersect(targets_in_ChIP,regulon_targets))  # OBSERVED SUCCESS

        # phyper(q=q_hyper, m=m_hyper,n=n_hyper,k=k_hyper,lower.tail = F, log.p = FALSE)

        # FISHER TEST: 
             # {{notation from wikipedia: Hypergeometric(N,K,n) with p(k)}} DIFFERENT FROM ABOVE!!
                # N_hyper <- length(total_universe_targets) # UNIVERSE

                # # targets_in_ChIP <- unique(ChipPeak_results[ChipPeak_results$peakName==tf,c("gene_name")]) # SUCCESS
                #  K_hyper <- length(targets_in_ChIP) # Number of SUCCESS
                 
                #  n_hyper <- length(regulon_targets) # number of DRAWNS

                #  k_hyper <-length(intersect(regulon_targets,targets_in_ChIP)) # Number of observed SUCCESS

        # For Teatesting:
        GRNinChip <- length(intersect(regulon_targets,targets_in_ChIP)) # regulon_results$target)) # goods
        
        GRNNoChip <- sum(! regulon_targets %in% targets_in_ChIP) # mistakes: n_hyper - k_hyper

        chipNoGRN   <- sum(!targets_in_ChIP %in% regulon_targets) # K_hyper - k_hyper

        NoChipNoGRN <- sum(! total_universe_targets %in% union(targets_in_ChIP,regulon_targets)) # regulon_results[[tf]] # Should be N_hyper + k_hyper - n_hyper - K_hyper # but does not match
        # NoChipNoGRN <- length(setdiff(total_universe_targets, union(targets_in_ChIP,regulon_targets)))
         
    TeaTasting <- matrix(c(GRNinChip, chipNoGRN , GRNNoChip , NoChipNoGRN ),
       nrow = 2,
       dimnames = list(GRN = c("Edge", "NoEdge"),
                       ChIPseq = c("Edge", "NoEdge")))

fisher_test <- fisher.test(TeaTasting,alternative=alternative)
# odds_res <- fmsb::oddsratio(GRNinChip, chipNoGRN, GRNNoChip, NoChipNoGRN, p.calc.by.independence = TRUE  )             
Results_phyper <- data.frame(   GRNinChip=GRNinChip,
                                GRNNoChip=GRNNoChip,
                                chipNoGRN=chipNoGRN,
                                NoChipNoGRN=NoChipNoGRN,
                                signif=ifelse(fisher_test$p.value<0.05,1,0),
                                p_val_Fisher_test=fisher_test$p.value,
                                conf.int_down=fisher_test$conf.int[1],
                                conf.int_up=fisher_test$conf.int[2],
                                odds_ratio=fisher_test$estimate,
                                row.names=tf)
if(out_targets){
    TF <- rep(tf,sum(GRNinChip,GRNNoChip))
    ChIPstatus_str <- c(rep("GRNinChip",GRNinChip),rep("GRNNoChip",GRNNoChip))
    ChIPstatus_bool <- c(rep(1,GRNinChip),rep(0,GRNNoChip))
    targets <- c(intersect(regulon_targets,targets_in_ChIP),setdiff(regulon_targets,targets_in_ChIP))
    targets_results_df=data.frame(targets=targets,ChIPstatus_str=ChIPstatus_str,ChIPstatus_bool=ChIPstatus_bool, TF=TF)
    return(list(Results_phyper=Results_phyper,
                targets_results_df=targets_results_df))
}

return(Results_phyper)
}

### LOAD RESULTS DATA-----------------------------------------------

#  PROMOTE/GRN DATA
print("Loading results")
targets_list <- rjson::fromJSON(file = '~/input_data/targets_list.json')
gene_info <- read.csv('~/input_data/promote_v1.gene_info.txt', sep="\t")
TFs <- gene_info$uniq_isos[gene_info$regulator==1]
print(paste("Total number of TFs in Human Protein Atlas:", length(TFs)))
promote_genes <- rownames(read.csv('~/input_data/mrm_sm_3500.txt', sep="\t"))
print(paste("TFs present in PROMOTE dataset:",length(intersect(TFs, promote_genes))))
total_universe_targets <- promote_genes[!promote_genes%in%TFs]
print(paste("Targets present in PROMOTE dataset (remaining):",length(total_universe_targets)))



regulon_TFs <- names(readRDS('~/output/sm_3500_50b_regulons.rds')$regulons)
regulon_results <- read.csv('~/output/tables/allTFs_diff_importance.csv') # 250k edges // 22k
regulon_results <- regulon_results[regulon_results$TF%in%regulon_TFs,] # 85k edges
# We will select the top 25% edges (for GRNboost2) from TF present in our rewired regulons




## See CHIPSEQ.R/CHIPSEQ_PARALLEL.r scripts to see how these objects were generated
all_targets_ChipPeak <- readRDS('~/output/all_PROMOTE_targets_intersect_ChipPeak.rds')
# validated_Tfs_ChipPeak <- readRDS('./validated_targets_ChipPeak.rds')


results_to_use <- all_targets_ChipPeak
ChipPeak_results <- results_to_use[results_to_use$gene_name%in%total_universe_targets,]

# TFs_intersect_rewired_regulons&ChipDataAnalyzed
TFs_intersect <- intersect(unique(regulon_results$TF),unique(results_to_use$peakName)) 

# get_stats_chipseq(tf,ChipPeak_results,regulon_results)

# out <- lapply(TFs_intersect,get_stats_chipseq,ChipPeak_results,regulon_results)
# chip_stats_results <- do.call(rbind,out)
# head(chip_stats_results)
# write.table(chip_stats_results,'./all_PROMOTE_targets_intersect_ChipPeak_stats.tsv',sep="\t")



out <- lapply(TFs_intersect,get_stats_chipseq,ChipPeak_results,regulon_results,"greater",T)
chip_stats_results <- lapply(out,function(x) x[[1]])
chip_stats_results <- do.call(rbind,chip_stats_results)
head(chip_stats_results)
write.table(chip_stats_results,'~/output/tables/GRNboost_regulon_top_diff_targets_intersect_ChipPeak_stats.tsv',sep="\t")
targets_with_Chip <- lapply(out,function(x) x[[2]])
targets_with_Chip <- do.call(rbind,targets_with_Chip)

print(paste(round(100* mean(chip_stats_results$signif),1),"% of rewired regulons have significant target enrichment on ChIPseq data"))
head(targets_with_Chip)


xlsx::write.xlsx(chip_stats_results,'~/output/tables/Supplementary_data_ChIPseq_GRNboost.xlsx',sheetName="Enrichment",row.names = T)
xlsx::write.xlsx(targets_with_Chip,'~/output/tables/Supplementary_data_ChIPseq_GRNboost.xlsx',sheetName="GRN rewired regulon with ChipSeq evidence",append =TRUE, row.names=F)
print("Saving results in: Supplementary_data_ChIPseq.xlsx")
print("DONE SCRIPT!")
# targets_with_Chip allows to guet most enriched TFs

# TeaTasting <- matrix(c(TraReinChip, chipNoTraRe , TraReNoChip , NoChipNoTraRe ),
#        nrow = 2,
#        dimnames = list(TraRe = c("Edge", "NoEdge"),
#                        ChIPseq = c("Edge", "NoEdge")))

