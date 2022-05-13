# SCRIPT TO EXTRACT THE SUPPLEMENTARY TABLES OF REGULON EXAMPLE
regulons <- readRDS( './output/sm_3500_50b_regulons.rds')

## NOTE: For the primer acquisition new names were used (Shown in tables and figures) 

# Table for name conversion: 
#  TF   New name    Old name
# ZNF3  "LAMTOR4"   "C7orf59"
# MXD1  "DMTN"      "EPB49" 
# MXD1  "MSRB1"     "SEPX1"
# MYB   "MSRB1"     "SEPX1"

# Targets names (old name) used in the qRT-PCR assay
targets_list <- list(
  ELK3= c("SEPP1","FAM65A","TAX1BP3","NCK1","SH3BP5","ARL15","ACTN4","ARL10","RAPGEF1","PIK3CA","MAP3K12","CLASP1","DUSP18","MYO18A","SRGAP2P2","RAB11B"),
  MXD1 = c("KEL","GCA","EPB49","SEPX1","MMP25","PTPN6","HBQ1"),
  MYB = c("GCA","ITGA2B","ERMAP","NRGN","SEPX1","FAM178B","HBQ1"),
  ZNF91=c("RSBN1","GEMIN5","MACROD2"),
  ZNF3 = c("C7orf59","LRRC27","DTX2","PSMG3","NUDCD3","STYXL1")
)


print_regulons <- function(targets_list,regulons){
  TFs <- names(targets_list)
  df <- list()
  for(TF in TFs){
    idx <- which(regulons$regulons[[TF]]$target%in%targets_list[[TF]])  
    df[[TF]] <- regulons$regulons[[TF]][idx,]
  }
  df_out <- do.call(rbind,df)
  rownames(df_out) <- NULL
  df_out$p_value <- format(df_out$p_value, scientific = T,digits = 3)
  return(df_out)
}

selected_regulons <- print_regulons(targets_list,regulons)

# All included!
table(selected_regulons$driver)%in%sapply(targets_list,length)
# Change names:
# Table for name conversion: 
#  TF   New name    Old name
# ZNF3  "LAMTOR4"   "C7orf59"
# MXD1  "DMTN"      "EPB49" 
# MXD1  "MSRB1"     "SEPX1"
# MYB   "MSRB1"     "SEPX1"

selected_regulons$target[selected_regulons$target=="C7orf59"] <- "LAMTOR4"
selected_regulons$target[selected_regulons$target=="EPB49"] <- "DMTN"
selected_regulons$target[selected_regulons$target=="SEPX1"] <- "MSRB1"

write.csv(selected_regulons,file="./output/tables/suppl_tables_S11-S15.csv", row.names = F)

################################ EOF #############################################
