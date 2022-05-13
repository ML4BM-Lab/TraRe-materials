#### SCRIPT TO CREATE SUPPLEMENTARY EXCEL FILE

# Read all the files ------------------------------------------------------

## Community enrichment GO

GO_BP_enrichment_regmods <- readRDS('./output/other/GO_BP_enrichment_regmods.rds')
regmod_go_enrichment <- GO_BP_enrichment_regmods@compareClusterResult

# TF GO ENRICHMENT

tf_go_enrichment <- read.csv('./output/tables/GO_enrichment_tfs.csv')


## REGULON RESULTS

regulon_object <- readRDS('./output/sm_3500_50b_regulons.rds')

regulon_results <- do.call(rbind,regulon_object$regulons)

## REGULON ENRICHMENT
GO_BP_cluster_regulons <- readRDS('./output/other/GO_enrichment_regulons_0.0001.rds')
regulon_go_enrichment <- GO_BP_cluster_regulons@compareClusterResult



# Export to excel ---------------------------------------------------------

final_list <- list(TableS6=regmod_go_enrichment,
                   TableS7=tf_go_enrichment,
                   TableS8=regulon_results,
                   TableS9=regulon_go_enrichment)

complete_excel <- openxlsx::createWorkbook()

sapply(names(final_list),function(x){
  openxlsx::addWorksheet(wb=complete_excel,sheetName=x)
  openxlsx::writeData(complete_excel,
                      sheet = x,
                      final_list[[x]],
                      startRow = 1,
                      startCol = 1,
                      rowNames = FALSE,
                      keepNA = TRUE)
})

openxlsx::saveWorkbook(complete_excel,file = './output/tables/supplementary_tables_S6-S9.xlsx')

################################ EOF #############################################