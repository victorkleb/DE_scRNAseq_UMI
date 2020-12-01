# method = c("edgeR", "DESeq2", "limma-trend", "limma-voom")

library(SummarizedExperiment)

library ( muscat )

limma_voom_DE_dsn <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/limma_voom_DE_B_cells.csv"

sce_dsn <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/sce.RData"
load ( sce_dsn)



pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)



# pseudobulks for 1st subpopulation
t(head(assay(pb)))



res_limma_voom <- pbDS(pb, verbose = TRUE, method = "limma-voom"  )   # explictly state method

# access results table for 1st comparison
tbl <- res_limma_voom$table[[1]]
# one data.frame per cluster
names(tbl)




# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))


write.csv( k1, file = limma_voom_DE_dsn )
