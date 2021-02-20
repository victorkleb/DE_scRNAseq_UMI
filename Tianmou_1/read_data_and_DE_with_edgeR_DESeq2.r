

#########################################################################################
#                                                                                       #    
#  program  read_data_and_DE_with_edgeR_DESeq2.r                                        # 
#                                                                                       #   
#########################################################################################

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

# data downloaded from https://github.com/Tianmou/scRNAseq-DE-comparison
# referenced in paper https://www.frontiersin.org/articles/10.3389/fgene.2019.01331/full



library( dplyr )
library ( edgeR )
library ( DESeq2 ) 

folder =  "D:/scRNA-seq_DE/simulation_Tianmou_1/"


# outputs 
 
counts_high_dsn   <-  paste0 ( folder, "counts_high.csv" ) 
counts_low_dsn    <-  paste0 ( folder, "counts_low.csv" ) 
trueDE_dsn        <-  paste0 ( folder, "trueDE.csv" ) 
 
DE_low_dsn <-  paste0 ( folder, "edgeR_and_DESeq_stats_low.csv" ) 
DE_high_dsn <-  paste0 ( folder, "edgeR_and_DESeq_stats_high.csv" ) 
 
 
# input - downloaded from Github
 
sim_data_dsn <-  paste0 ( folder, "BPsimData_0.05DE_lfc1.Rdata" )

#########################################################################################

load (sim_data_dsn)

print (class(data))
print (dim(data))
print (data[1:5, 1:5])


# input results are real, convert to integer
data = floor(data) 
storage.mode(data) <- "integer"
print (data[1:5, 1:5])

print (length( iso.low ))
print (head ( iso.low ))

print (length ( iso.high ))
print (head ( iso.high ))

print (length( trueDE ))
print (head ( trueDE ))


data_low = data [ iso.low, ]
print (dim(data_low))
print ( data_low[1:5, 1:5])

data_high = data [ iso.high, ]
print (dim(data_high))
print ( data_high[1:5, 1:5])


write.csv( data_low, file = counts_low_dsn )
write.csv( data_high, file = counts_high_dsn )
write.csv( trueDE, file = trueDE_dsn )

#########################################################################################

# did not this find stated in paper, but EDA suggests this is the correct split between cell clusters

cell_groups = rep(1:2, each=80 ) 
group_data = factor ( cell_groups )
print ( group_data )



df_edgeR  = edgeR ( data_low, group_data )
df_DESeq2 = DESeq2  ( data_low, group_data )

df_DE_low = merge ( df_edgeR, df_DESeq2, by=0, all=TRUE)
names(df_DE_low)[1]<- 'Gene' 
print ( head( df_DE_low ) )
write.csv( df_DE_low, file = DE_low_dsn )



df_edgeR  = edgeR ( data_high, group_data )
df_DESeq2 = DESeq2  ( data_high, group_data )

df_DE_high = merge ( df_edgeR, df_DESeq2, by=0, all=TRUE)
names(df_DE_high)[1]<- 'Gene' 
print ( head( df_DE_high ) )
write.csv( df_DE_high, file = DE_high_dsn )

