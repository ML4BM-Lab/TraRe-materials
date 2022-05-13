# Script to plot Figure 3B of enrichment of rewired TFs (TraRe's manuscript)

# REQUIRED LIBRARIES LOADINGS ----------------------------------------------------------------
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# versions: 
# clusterProfiler_4.0.5
# org.Hs.eg.db_3.13.0 
# AnnotationDbi_1.54.1
# DOSE_3.18.3   
# ggplot2_3.3.5 
# forcats_0.5.1
# stringr_1.4.0 
# dplyr_1.0.7

# Automatically attached:
# IRanges_2.26.0        
# S4Vectors_0.30.2
# Biobase_2.52.0       
# BiocGenerics_0.38.0 

# load file
rewired_drivers <- rownames(read.delim2('./output/tables/table_rewired_tfs.tsv'))
# ENRICHMENT --------------------------------------------------------------
# Transform labels
# Annotation
OrgHs <- org.Hs.eg.db
# Gene nomenclature
gene_entrez <- AnnotationDbi::mapIds(OrgHs,
                      keys= rewired_drivers,
                      column="ENSEMBL",
                      keytype="SYMBOL",
                      multiVals="first")

# enrichGO: GO over-representation test (Hypergeometric)
biol <- clusterProfiler::enrichGO(gene = gene_entrez,
                 keyType = "ENSEMBL",
                 OrgDb = OrgHs,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

biol <- mutate(biol, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
biol <- mutate(biol, FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))

# Save results
write.csv(biol@result,'./output/tables/GO_enrichment_tfs.csv', row.names = F )

# With ggplot2
df <- biol@result
df <- df[df$p.adjust<0.01,]

pdf('./output/figures/GO_TFs_dotplot.pdf', width=7)
ggplot(df, showCategory=20,aes(x = Count,y = forcats::fct_reorder(Description, Count))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = DOSE::parse_ratio(GeneRatio))) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))+
  scale_x_discrete(limits=c(2,4,6,8,10)) +
  theme_bw() + 
  labs(size = "GeneRatio")+
  xlab("Count") +
  ylab(NULL) + 
  ggtitle("Biological Process Ontology enrichment")
dev.off()
################################ EOF #############################################