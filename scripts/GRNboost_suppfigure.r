#!/usr/bin/env Rscript
#' ######################################################
#' ######################################################
#'              GNRboost2 results analysis for
#'          Cancer Research TraRe manuscript revision
#' ######################################################
#' ######################################################
library(ggplot2)
library(dplyr)

targets_list_path= './GRNboost2/targets_list.json'
targets_list <- rjson::fromJSON(file = targets_list_path)

GRN_path_R = './GRNboost2/network_responders.tsv'
grnboostR <- read.csv(GRN_path_R, sep = ",")[, -1]

GRN_path_NR = './GRNboost2/network_non_responders.tsv'
grnboostNR <- read.csv(GRN_path_NR, sep = ",")[, -1]

matched_edge_R <- read.csv('./GRNboost2/R/validating_edges_importance_pvals.csv')
matched_edge_NR <- read.csv('./GRNboost2/NR/validating_edges_importance_pvals.csv')


matched_edge_df <- rbind(matched_edge_R,matched_edge_NR)
matched_edge_df$pheno <- factor(c(rep("R",nrow(matched_edge_R)),rep("NR",nrow(matched_edge_NR))),levels=c("R","NR"))

edge_importancesR <- grnboostR[which(grnboostR$TF %in% names(targets_list)),]
edge_importancesNR <- grnboostNR[which(grnboostNR$TF %in% names(targets_list)),]
pheno <- factor(c(rep("R",nrow(edge_importancesR)),rep("NR",nrow(edge_importancesNR))), levels=c("R","NR"))

edges_all <- rbind(edge_importancesR,edge_importancesNR)
edges_all <- cbind(edges_all, pheno)
final_df <- left_join(edges_all, matched_edge_df, by=c("TF"="TF","target"="target","importance"="importance","pheno"="pheno") )


values_color <- c("#4284c1","#c08541")
histograms <- ggplot(final_df, aes(x=importance,y=..ncount..))+
geom_histogram(aes(fill=pheno))+
geom_vline(data=filter(matched_edge_df, sig==1),aes(xintercept=importance), colour="red")+
geom_vline(data=filter(matched_edge_df, sig==0),aes(xintercept=importance),linetype="dotdash", colour="black")+
scale_x_log10()+
scale_fill_manual(values=values_color)+
guides(fill=guide_legend(title="Response"))+
facet_grid(rows = vars(TF), cols=vars(pheno),scales = 'free_x')+
theme_bw()+
theme(text = element_text(size = 10))

diff_edges_R <- setdiff(edge_importancesR[,-3], edge_importancesNR[,-3])
diff_edges_NR <- setdiff(edge_importancesNR[,-3], edge_importancesR[,-3])
non_found <- sum(nrow(diff_edges_R),nrow(diff_edges_NR))

matched_edges <- full_join(grnboostR, grnboostNR,
                        by=c("TF"="TF","target"="target"),
                        suffix = c("_R", "_NR") )%>% 
                mutate(diff_importance=abs(importance_R-importance_NR))

matched_edges$diff_importance[is.na(matched_edges$diff_importance)] <- 0

matched_edges <- matched_edges[order(matched_edges$diff_importance,decreasing = TRUE),]

final_matched_edges <- matched_edges[matched_edges$diff_importance!=0,]

ven <- lapply(names(targets_list), function(driver){
ggvenn::ggvenn(list(R=edge_importancesR[edge_importancesR$TF==driver,2],
                    NR=edge_importancesNR[edge_importancesNR$TF==driver,2]),
            fill_color = values_color,   
            stroke_size = 0.5, set_name_size = 4,text_size=1.7)

})
vens <- cowplot::plot_grid(ven[[1]],ven[[2]],ven[[3]],ven[[4]],ven[[5]],ncol = 1)
pdf('./GRNboost2/venns_grid.pdf', width=2,height = 10)
vens
dev.off()

pdf('./GRNboost2/histogram_grid.pdf', width=7,height = 10)
histograms
dev.off()

# outR <- sapply(names(targets_list), get_driver_distribution, matched_edge[[1]], grnboostR, log=TRUE)
# outNR <- sapply(names(targets_list), get_driver_distribution, matched_edge[[2]], grnboostNR, log=TRUE)
print("DONE SCRIPT!")
