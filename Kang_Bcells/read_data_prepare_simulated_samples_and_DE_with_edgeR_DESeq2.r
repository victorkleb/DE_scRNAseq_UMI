################################################################################################################################ 
##                                                                                                                            ##
##  start program: read_data_prepare_simulated_samples_and_DE_with_edgeR_DESeq2.r                                             ##
##                                                                                                                            ## 
################################################################################################################################

edgeR  <- function ( counts_array, group_data )
{
  DGE_List = DGEList ( counts=counts_array, group=group_data )
  DGE_List = calcNormFactors(DGE_List)

  design_anova = model.matrix( ~group_data, data=DGE_List$samples ) 

  DGE_List_power = estimateGLMTrendedDisp  (DGE_List, design_anova, method="power" )
  fit = glmQLFit( DGE_List_power, design_anova )
  qlf = glmQLFTest( fit, coef=2:2 )
  edgeR_results = qlf$table

  FDR <- p.adjust(qlf$table$PValue, method="BH")
  edgeR_results$edgeR_FDR = FDR
  
  df_edgeR = dplyr::select ( edgeR_results, F, PValue, edgeR_FDR )
  names ( df_edgeR )[1] = "edgeR_F"
  names ( df_edgeR )[2] = "edgeR_PValue"  
  
  return ( df_edgeR )
}  



DESeq2 <- function ( counts_array, group_data )
{
  counts_array_p1 = counts_array + 1

  dds <- DESeqDataSetFromMatrix(countData = counts_array_p1,
                                  DataFrame ( group_data ),
                                  design= ~ group_data )
								
  dds <- DESeq(dds )
  DESeq2_results <- results(dds)
  df_DESeq2_results = as.data.frame( DESeq2_results )
 
  df_DESeq2 = dplyr::select ( df_DESeq2_results, stat, padj )
  names ( df_DESeq2 )[1] = "DESeq2_stat"
  names ( df_DESeq2 )[2] = "DESeq2_padj" 
  
  return ( df_DESeq2 )
}  

#########################################################################################

# http://www.bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html

# follow vignette to obtain B cell counts for highly expressed genes


library ( ExperimentHub )
library ( scater )
library ( splatter )
library ( edgeR )
library ( DESeq2 )




number_of_simulations = 10

folder =  "D:/scRNA-seq_DE/simulation_Kang_Bcells/"

# output base dsns 

splatter_counts_base   <-  paste0 ( folder, "splatter_counts_" ) 
splatter_rows_base     <-  paste0 ( folder, "splatter_rows_" )
splatter_columns_base  <-  paste0 ( folder, "splatter_columns_" ) 
splatter_DE_base      <-  paste0 ( folder, "splatter_DE_" ) 

################################################################################################################################ 
 
eh <- ExperimentHub()
query(eh, "Kang")


(sce <- eh[["EH2259"]])

print ( dim(sce) )

sce <- sce[rowSums(counts(sce) > 0) > 0, ]
print ( dim(sce) )


qc <- perCellQCMetrics(sce)


ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
print ( dim(sce) )


sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
print ( dim(sce) )


# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

################################################################################################################################

### select highly expressed B cells

B_cells = sce [ , sce$cell=='B cells' ]
print ( dim( B_cells ) )
rm ( sce  )


# 500 genes expressed in the most cells
cells_expressed = rowSums(counts(B_cells) > 1)
cells_expressed_sorted = sort ( cells_expressed, decreasing=TRUE )
print ( cells_expressed_sorted[500] )

loc_select = ( cells_expressed>= cells_expressed_sorted[500] )
B_cells_select = B_cells[ cells_expressed>= cells_expressed_sorted[500] ]
print ( dim ( B_cells_select ) )
# and check
print ( min ( rowSums(counts(B_cells_select) > 1) ) )
rm ( B_cells  )


################################################################################################################################

# https://bioconductor.org/packages/release/bioc/vignettes/splatter/inst/doc/splatter.html

# splatter: compute parameters from genuine data

count_data = counts( B_cells_select )
counts_array = data.matrix ( count_data )

# Check that counts_array is an integer matrix
print ( class(counts_array) )
print ( typeof(counts_array) )

# Check the dimensions, each row is a gene, each column is a cell
print ( dim(counts_array) )

# Show the first few entries
counts_array[1:5, 1:5]


# splatter parameters 
params0 <- splatEstimate(counts_array)

 


# splatter simulations and DE with edgeR and DESeq2: results to csv files

for ( simulation in 1:number_of_simulations )
{  
  str_sim = toString ( simulation )
  print ( paste ( "simulation: ", str_sim ) )

  params <- setParams(params0, seed=10000 + simulation )
  sim <- splatSimulateGroups(params, group.prob = c(0.5, 0.5), de.prob = 0.3 )


  col_data = colData( sim )
  row_data = rowData( sim )
  simulated_count_data = counts( sim )
  count_matrix = as.matrix( simulated_count_data )


  df_edgeR  = edgeR   ( count_matrix, col_data$Group )
  df_DESeq2 = DESeq2  ( count_matrix, col_data$Group )

  df_DE = merge ( df_edgeR, df_DESeq2, by=0, all=TRUE)
  names(df_DE)[1]<- 'Gene' 
  print ( head( df_DE ) )
  



  row_data_dsn    <- paste0 ( splatter_rows_base,   str_sim, ".csv" )
  column_data_dsn <- paste0 ( splatter_columns_base, str_sim, ".csv" )
  counts_out_dsn  <- paste0 ( splatter_counts_base, str_sim, ".csv" )
  DE_dsn          <-  paste0 ( folder, "splatter_DE_", str_sim, ".csv" )
  
  write.csv( col_data, file = column_data_dsn )
  write.csv( row_data, file = row_data_dsn )
  write.csv( count_matrix, file = counts_out_dsn )
  write.csv( df_DE, file = DE_dsn )
}

 