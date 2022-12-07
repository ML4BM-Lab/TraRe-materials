
#'##############################################################################
#'#########################	 ChIPseq analysis with ChIPpeakAnno ################
#'############################## IN PARALLEL COMPUTATION #######################
#' Tutorial: https://www.bioconductor.org/packages/devel/bioc/
#' 			vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html
#' 
#' We use a Public DB for ChIPseq data Remap2022. It contains multiple data
#' 			Have a look at (https://remap.univ-amu.fr/download_page).
#' We used the 'Non redundant peaks' (counted peaks for each TF) for all
#' the experiments in Homo sapiens (mix of tissues, experiments).
#' 
#' We use annotatePeakInBatch to annotate the Remap peaks with the genomic
#' features in the AnnotationData (genecode.v33 in this case) within certain distance
#' away specified by maxgapor and other parameters.
#' SCENARIOS for match:
#' a) if  the region is upstream the gene start, no further than 5kb.
#' b) if the region maps into some internal regions like introns but not too
#' 		distant to the nearest feature, a maximum of 5kb.
#' c) if the region overlaps with the start of the described gene or feature.
#' d) if the region is a smaller region of exon or described region of the gene
#' 
#' We use all the remap peaks for a specific TF (not only our predicted target-regions)
#' to be able to compute enrichment tests.
#' 
#' Because we want this information for each TF we parallelize the computation.
#' 
#' #############################################################################
#'  Author: Irene Marin (based on Guillermo Serrano's script)
#' #############################################################################

### ------LOAD DATA---------------------------------
## See CHIPSEQ.R script to see how these objects were generated

print("Loading data")
gene_info <- read.csv('~/input_data/promote_v1.gene_info.txt', sep="\t")
TFs <- gene_info$uniq_isos[gene_info$regulator==1]

# all_peaks <- setNames(read.table('remap2022_nr_macs2_hg38_v1_0_TFs.bed', stringsAsFactors = F), c('chr', 'start', 'end','peakName', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'))
# saveRDS(all_peaks, "./remap2022_chippeaks.rds")

all_peaks <- readRDS("~/input_data/remap2022_chippeaks.rds")
all_TFs<- unique(all_peaks$peakName)

TFs_intersect <- intersect(TFs,all_TFs)
length(TFs_intersect)

gencode33all <- readRDS('~/input_data/DF_gencode.v33.annotation_filtered.rds')
G19AnnotChIP_RL <- readRDS('~/input_data/GRanges.gencode.v33.annotation_filtered.rds')

# Set number of cores

NrCores <- parallel::detectCores(logical = FALSE)/2
parallClass <- BiocParallel::bpparam()
parallClass$workers <- NrCores


FuncToParallelize <- function(tf, peakName){
	result_idx <- which(grepl(pattern=paste0("^",tf,"$"), peakName))
}
print("First parallelization")
result1_idx <- BiocParallel::bptry({ BiocParallel::bplapply(TFs_intersect, FuncToParallelize, all_peaks$peakName, BPPARAM = parallClass)})
names(result1_idx) <- TFs_intersect
result_idx_to_check <- result1_idx[which(sapply(result1_idx,length)!=0)]
length(result_idx_to_check)


FuncToParallelize2 <- function(result_idx_to_check, all_peaks, G19AnnotChIP_RL,gencode33all){
	
	peaks_tmp <- all_peaks[result_idx_to_check,]
	peaks_tmp$strand <- '+'

	peaks_bed <- ChIPpeakAnno::toGRanges(peaks_tmp, format="BED", header=TRUE)
	names(peaks_bed) <- seq(1, peaks_bed@elementMetadata@nrows)

	beddta_Annot_df<- as.data.frame(ChIPpeakAnno::annotatePeakInBatch(peaks_bed, AnnotationData = G19AnnotChIP_RL, output = "both", multiple = F, maxgap = 0))

	beddta_Annot_df_f5 <- beddta_Annot_df[((beddta_Annot_df$fromOverlappingOrNearest == "NearestLocation") &
											(beddta_Annot_df$insideFeature == "upstream") & 
											(abs(beddta_Annot_df$distancetoFeature)<5000)) | # scenario a)
											((beddta_Annot_df$fromOverlappingOrNearest == "NearestLocation") &
											(beddta_Annot_df$insideFeature == "inside") & 
											(abs(beddta_Annot_df$distancetoFeature)<5000)) | # scenario b)
											((beddta_Annot_df$fromOverlappingOrNearest == "NearestLocation") & 
											(beddta_Annot_df$insideFeature == "overlapStart")) | # scenario c)
											((beddta_Annot_df$fromOverlappingOrNearest == "NearestLocation") & 
											(beddta_Annot_df$insideFeature == "includeFeature")),] # scenario d)
	
	# Add the gene_name (symbol ID) using the Annotation data
    beddta_Annot_df_f5 <- merge(beddta_Annot_df_f5, gencode33all[,c('gene_id', 'gene_name')], by.x='feature', by.y='gene_id')
	return(beddta_Annot_df_f5)
}
print("Second parallelization")
result <- BiocParallel::bptry({ BiocParallel::bplapply(result_idx_to_check, FuncToParallelize2, all_peaks, G19AnnotChIP_RL, gencode33all, BPPARAM = parallClass)})
final_result <- do.call(rbind,result)


print(paste("Total TFs with no peaks", sum(sapply(result1_idx,length)==0)))


saveRDS(final_result,'~/output/all_PROMOTE_targets_intersect_ChipPeak.rds')
write.table(final_result,'~/output/all_PROMOTE_targets_intersect_ChipPeak.tsv', sep="\t")
print(paste("Results saved in rds and tsv as all_PROMOTE_targets_intersect_ChipPeak! "))
print("DONE SCRIPT!")