 

library ( glmpca )
library ( mclust )


glmpca_results_dsn = "D:/scRNA-seq_DE/Zhengmix8eq/RNA seq data - glmpca 10 output factors poisson - 195 genes.csv"
Mclust_out_dsn = "D:/scRNA-seq_DE/Zhengmix8eq/RNA seq data - Mclust output - 2-20 clusters from glmpca 10 output factors poisson - 195 genes.csv"

count_dsn    <-  "D:/scRNA-seq_DE/Zhengmix8eq/filtered_UMI_counts.csv"
count_data <- as.matrix ( read.csv ( count_dsn, header=TRUE, row.names="X"  ) )


################################################################################################
			
glmpca_res = glmpca ( count_data, 10, fam="poi"  )
factors = glmpca_res$factors
sample_ID = row.names( factors)
rownames(factors) <- c()



df_clusters = data.frame ( sample_ID )

for  ( num_clusters in 2:20 )
{		
Mclust_res = Mclust( factors, num_clusters,verbose=TRUE)
clusters = Mclust_res$classification
str_num_clusters = toString (num_clusters)
df_clusters$add = clusters 
colnames(df_clusters)[num_clusters] = str_num_clusters
}



write.csv( factors, file = glmpca_results_dsn )			
write.csv( df_clusters, file = Mclust_out_dsn )																	 
								 