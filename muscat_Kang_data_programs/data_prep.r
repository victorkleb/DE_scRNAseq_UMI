################################################################################################################################ 
##                                                                                                                            ##
##  start program: data_prep.r                                                                                                ##
##                                                                                                                            ## 
################################################################################################################################

# http://www.bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html
# https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

library(ExperimentHub)
library (muscat)



#counts_dsn      <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/sce_counts.csv"
#row_data_dsn    <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/sce_row_data.csv"
#column_data_dsn <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/sce_column_data.csv"

sce_dsn <-  "D:/scRNA-seq_DE/Crowell_Kang18_8vs8/sce.RData"




eh <- ExperimentHub()
query(eh, "Kang")


(sce <- eh[["EH2259"]])

dim(sce)

sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)


ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)


sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)


sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce, 
    kid = "cell", # subpopulation assignments
    gid = "stim",  # group IDs (ctrl/stim)
    sid = "id",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids	
	
	
	
	
# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))	
	
	
# compute UMAP using 1st 20 PCs
sce <- runUMAP(sce, pca = 20)
	
	
# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_minimal() + theme(aspect.ratio = 1)	
	
# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (dr in c("TSNE", "UMAP"))
  for (col in c("cluster_id", "group_id"))
    .plot_dr(sce[, cs100], dr, col)	
	
	




#col_data = colData( sce )
#row_data = rowData( sce )
#count_data = counts( sce )
#count_matrix = as.matrix( count_data )


#write.csv( col_data, file = column_data_dsn )
#write.csv( row_data, file = row_data_dsn )
#write.csv( count_matrix, file = counts_dsn )



save( sce, file = sce_dsn)



################################################################################################################################ 
##                                                                                                                            ##
##  end program: data_prep.r                                                                                                  ##
##                                                                                                                            ## 
################################################################################################################################