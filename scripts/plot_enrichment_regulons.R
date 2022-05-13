# Script to plot Figure 3D of enrichment of regulons (TraRe's manuscript)

# REQUIRED LIBRARIES LOADINGS ----------------------------------------------------------------
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

suppressPackageStartupMessages(library(rrvgo))
suppressPackageStartupMessages(library(GOSemSim))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(hues))

# versions: 
# clusterProfiler_4.0.5
# org.Hs.eg.db_3.13.0 
# AnnotationDbi_1.54.1 
# rrvgo_1.4.4 
# GOSemSim_2.18.1   
# ggplot2_3.3.5 
# forcats_0.5.1
# hues_0.2.0 

# Automatically attached:
# IRanges_2.26.0
# S4Vectors_0.30.2
# Biobase_2.52.0       
# BiocGenerics_0.38.0      

# OBJECTS
regulons_p <- './output/sm_3500_50b_regulons.rds'
regulons <- readRDS(regulons_p)

# ENRICHMENT --------------------------------------------------------------
# Transform labels
# Annotation
OrgHs <- org.Hs.eg.db
genelist <- lapply(regulons$regulons, function(x) AnnotationDbi::mapIds(OrgHs, column="ENTREZID",
                                                keytype="SYMBOL",
                                                multiVals="first",
						                                    keys=c(x$driver[1],x$target)))

# Compare cluster biological process
GO_BP_cluster_regulons <- clusterProfiler::compareCluster(geneCluster = genelist,
                     fun = clusterProfiler::enrichGO,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     OrgDb=OrgHs,
                     pvalueCutoff = 0.0001,
                     qvalueCutoff = 0.0001,
                     maxGSSize = 2000, 
                     readable=TRUE)

# Save
saveRDS(GO_BP_cluster_regulons,'./output/other/GO_enrichment_regulons_0.001.rds')

# ENRICHMENT RESULT OBJECT
GO_BP_cluster_regulons <- readRDS('./output/other/GO_enrichment_regulons_0.001.rds')

# Reorder factor to show together TFs from same regulatory module
tfs_regulons <- names(regulons$regulons)

tfs_Coms <- c("CEBPE", "GATA1" , "KLF1" ,"MXD1" , "NFE2" , "FOXN1" ,"GLI1" , "MYB",
"ZNF91" ,"ZFHX4",
"PAX6", "SOX8",
"SNAI2" ,
"NR1H4",
"PRRX1" ,"KLF10", "MSC",
"ELK3" ,"DNMT1",
"ZNF174",
"RUNX1")

mix_tfs <- c(tfs_Coms,tfs_regulons[!tfs_regulons%in%tfs_Coms])

# Extract enrichment data.frame
df <- GO_BP_cluster_regulons@compareClusterResult

# Reorder clusters on columns (regulons) as factor 
mylevels <- mix_tfs[mix_tfs%in%unique(df$Cluster)]
df$Cluster <- as.character(df$Cluster)
df$Cluster <- forcats::fct_relevel(df$Cluster, mylevels)


# Calculate similar terms
d <- GOSemSim::godata(OrgHs, ont="BP")

simMatrix_all <- rrvgo::calculateSimMatrix(df$ID,
                                       orgdb = OrgHs,
                                       ont="BP",
                                       method="Rel",
                                       semdata = d)

reducedTerms_all <-  rrvgo::reduceSimMatrix(simMatrix_all,
                                threshold=0.75,
                                orgdb=OrgHs) # scores not used

# Add parent term column and reorder by it
idx <- match(df$ID,reducedTerms_all$go)
df$parentTerm <- reducedTerms_all$parentTerm[idx]
df_ordered <- df[order(df$parentTerm),]
df_ordered$Description <- factor(df_ordered$Description,levels=unique(df_ordered$Description))

# Some info messages
sprintf("Number of parent terms %d",dim(df_ordered)[1])
sprintf("Number of GO terms %d",length(unique(reducedTerms_all$parentTerm)))

# Prepare colors for plot
nterms <- length(unique(reducedTerms_all$parentTerm))
mycolors<- hues::iwanthue(n=nterms,hmin=1, hmax=350,cmin = 2,cmax = 180,lmin = 20,lmax=100)
names(mycolors) <- unique(df_ordered$parentTerm[order(df_ordered$parentTerm)])


# PLOT --------------------------------------------------------------------

pdf("./output/figures/Regulon_GO_clusterEnrichment_0.0001_thr0.75.pdf",width=15,height=20)
ggplot(df_ordered)+
geom_tile(aes(x=Cluster, y=Description ,fill=parentTerm, alpha=-log(p.adjust)))+
scale_y_discrete(limits = rev(levels(df_ordered$Description)))+
scale_fill_manual(values=mycolors)+
scale_alpha(range=c(0.5,1))+ # for clarity
theme_classic()+
theme(legend.direction="vertical",axis.text.y= element_text(size=5))+
guides(fill=guide_legend(ncol=1))
dev.off()

#######################################################################################

# To facilitate plot modifications with Adobe Illustrator here is the code to generate a similar 
# plot with pheatmap package to anotate columns and rows of the previous plot
# matrix for heatmap

myrows <- unique(df_ordered$Description)
mycols <- mylevels
mymatrix <- matrix(data=NA, nrow=length(myrows), ncol=length(mycols))
rownames(mymatrix) <- myrows
colnames(mymatrix) <- mycols
for (row in myrows ){
  for(col in mycols){
    x <- df_ordered[df_ordered$Cluster==col,]
    if(row%in%x$Description){
      mymatrix[row,col] <-  -log(df_ordered[df_ordered$Cluster==col & df_ordered$Description==row,"p.adjust"])
    }
  }
}

annotation_row <- c()
i=0
for (row in myrows){
  i <- i+1
  print(row)
  print(i)
  annotation_row <- c(annotation_row,unique( df_ordered[df_ordered$Description==row, "parentTerm"]))
}
annotation_row <- as.data.frame(annotation_row)
rownames(annotation_row) <- rownames(mymatrix)
colnames(annotation_row) <- "parentTerm"


annotation_col <- as.data.frame(c("M1","M1","M1","M1","M1",
                                  "M4","M5", "M6","M6","M7","M8","M9", "NC", "NC", "NC", "NC"))

rownames(annotation_col) <- as.character(mycols)
colnames(annotation_col) <- "Module"

mycolors<- hues::iwanthue(n=nterms,hmin=1, hmax=350,cmin = 2,cmax = 180,lmin = 20,lmax=100)
names(mycolors) <- unique(df$parentTerm[order(df$parentTerm)])

Module <- hues::iwanthue(n=length(unique(annotation_col$Module)))
names(Module) <- unique(annotation_col[,1])

col_list <- list(parentTerm=mycolors,Module=Module)

pdf("./output/others/Regulon_GO_clusterEnrichment_0.0001_thr0.75_heatmap.pdf",width=15,height=20)

pheatmap::pheatmap(mymatrix,na_col="white",
                   cluster_cols=F,cluster_rows=F,
                   fontsize = 5,
                   cellheigth=5, cellwidth = 10,
                   show_rownames = T,
                   annotation_row= annotation_row, 
                   annotation_colors = col_list,
                   annotation_col = annotation_col,
                   color = hcl.colors(50, "BluYl"))
dev.off()

################################ EOF #############################################