#!/usr/bin/env Rscript
#' ######################################################
#' ######################################################
#'              ARACNE-AP results analysis for
#'          Cancer Research TraRe manuscript revision
#'  Print histograms
#'  Create tables with edge pvalues
#' ######################################################
#' ######################################################

rm(list = ls())
# Get passed args
args = commandArgs(trailingOnly=TRUE)
GRN_path <- args[1]
targets_list_path <- args[2]
output_path <- args[3]

# Create dir if not exist already
if (file.exists(output_path)){
    output_path <- paste0(getwd(),'/',output_path)
} else {
    dir.create(paste0(getwd(),'/',output_path))
    output_path <- paste0(getwd(),'/',output_path)
}

print(paste("Saving output in",output_path))


print(paste("Import results from", GRN_path))
grnboost <- read.csv(GRN_path, sep = "\t")[,c(1,2,5)]
colnames(grnboost) <- c("TF","target","importance")
grnboost <- grnboost[order(grnboost$importance, decreasing = T),]
rownames(grnboost) <- NULL

print(paste("Import targets list from", targets_list_path,"\\n"))
targets_list <- rjson::fromJSON(file = targets_list_path)

### ------------------------------------------------------------
# Checks on GRN output if the validating edge exists
check_validated_edge <- function(driver, target, GRN) {
    driver_GRN <- GRN[GRN$TF == driver, ] # 
    rownames(driver_GRN) <- seq(1, nrow(driver_GRN)) 
    matched_edge <- driver_GRN[driver_GRN$target == target, ]
   if (nrow(matched_edge) == 0) {
        print(paste("No edge found for", driver, "and", target))
        return(data.frame(TF=driver,target=target,importance=NA))
    }
    return(matched_edge)
}

# Prints histogram of edges importances and the foudn matching validating edges
get_driver_distribution <- function(driver, matched_edge, GRN, log = F) { # nolint
    edge_importances <- GRN[which(GRN$TF == driver), 3]
    matched_importance <- matched_edge[[driver]]$importance
    info <- paste0(as.character(sum(!is.na(matched_importance))),"/",as.character(length(matched_importance)))

    proportion <- mean(!is.na(matched_importance))
    breaks <- length(edge_importances) /1.5
    if (log) {
        edge_importances <- log10(edge_importances)
        matched_importance <- log10(matched_importance)
        breaks <- length(edge_importances) / 100
    }


    hist(edge_importances,
        breaks = breaks,
        freq = T,
        xlab = "",
        main = paste("Edge importance distribution for", driver),
        sub = paste("Number of validating edges found:", info, "\n", "Total TF-target edges:", length(edge_importances))
    )
    abline(v = matched_importance, col = "red")
}

matched_edge <- list()
for (driver in names(targets_list)) {
    if(driver%in%grnboost$TF){
    matched_edge[[driver]] <- lapply(targets_list[[driver]], function(target) {
        check_validated_edge(driver, target, grnboost)
    })
    matched_edge[[driver]] <- do.call(rbind, matched_edge[[driver]])
    # rownames(matched_edge[[driver]]) <- NULL
    }
}


# get_driver_distribution(driver, matched_edge, grnboost)

pdf(paste0(output_path,'/validating_edges_histograms_fulldata.pdf'))
# get_driver_distribution(driver,matched_edge,grnboost)
 out <- sapply(names(targets_list) , get_driver_distribution, matched_edge, grnboost, log=TRUE)
dev.off()

### ------------------------------------------------------------
get_pvals_validating_targets <- function(driver=NULL, target=NULL, GRN=NULL,driver_matched_edge=NULL,alpha=0.05 ) {
    if (is.null(driver)){stop('driver missing')}
    if (is.null(target)){stop('target missing')}
    if (is.null(GRN)){stop('GRN missing')}
    if (is.null(driver_matched_edge)){stop('list_matched_edge missing')}
    if (alpha>0.05){warning('alpha is greater than 0.05!')}

    edge_importance <- driver_matched_edge[driver_matched_edge$target==target,3,drop=F]

    if (is.na(edge_importance)) {
        df <- data.frame(driver,target,NA,NA,NA, NA)
        colnames(df) <- c("TF","target","importance","edge_pval","sig","sig_pval")
    }else{
    edge_importances <- GRN[which(GRN$TF == driver), 3]
    sig_pval <- length(edge_importances)*alpha
    edge_pval  <-  as.numeric(rownames(edge_importance)) / length(edge_importances)
    # sig2 <- ifelse(edge_pval<alpha,1,0)
    sig <- ifelse(as.numeric(rownames(edge_importance))<sig_pval,1,0)
    df <- cbind(driver_matched_edge[driver_matched_edge$target==target,,drop=F],edge_pval,sig,sig_pval)
    }
    return(df)
}

# get_pvals_validating_targets(driver, target, grnboost)

pval_df <- do.call(rbind, lapply(names(targets_list), function(driver){
    out <- do.call(rbind, lapply(targets_list[[driver]], function(target) {
        get_pvals_validating_targets(driver,target,grnboost,driver_matched_edge=matched_edge[[driver]])}))
    out_ordered <- out[order(out$edge_pval),]
    return(out_ordered)}))

write.table(pval_df,paste0(output_path,'/validating_edges_importance_pvals.csv'),sep=",",row.names = F)

print("DONE!")
