


#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  MLE_cell_cluster_models.py                                            # 
#                                                                                       #   
#########################################################################################
#
#  invoke FUNCTION   MLE_cell_cluster_models   
#
#########################################################################################

import pandas as pd
import numpy  as np



import pickle
 
from pathlib import Path
import os



import sys
sys.path.append("D:/scRNA-seq_DE/DE_scRNAseq")


from  utilities import *
from FUNCTIONS_DE_scRNAseq  import *


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

pctl_list =  [ .1, .25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ]

#########################################################################################  


logfile_txt =  "MLE_cell_cluster_models.txt"

dict_MLE_results_pkl = "MLE_cell_clusters_all_genes.pkl"
NLL_LR_stat_summary_pkl = "NLL_LR_stat_summary_all_genes.pkl"
Mclust_clusterings_pkl = "Mclust_clusterings_all_genes.pkl"


UMI_counts_csv = "sce_full_Zhengmix8eq_counts.csv"
Mclust_clusterings_csv = "RNA seq data - Mclust output - 2-20 clusters from glmpca 10 output factors poisson - 195 genes.csv"


data_folder = Path ( r"D:/scRNA-seq_DE/Zhengmix8eq" )
data_path = Path ( data_folder )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle outputs
dict_MLE_results_dsn = data_path / dict_MLE_results_pkl
NLL_LR_stat_summary_dsn = data_path / NLL_LR_stat_summary_pkl
Mclust_clusterings_pkl_dsn = data_path / Mclust_clusterings_pkl


# CSV inputs
UMI_counts_dsn = data_path /  UMI_counts_csv
Mclust_clusterings_csv_dsn = data_path / Mclust_clusterings_csv

######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_dsn , index_col=0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )


##########  !!!  clusters from Mclust are numbered 1...N  not 0...N-1  
df_cell_clusters_in = pd.read_csv ( Mclust_clusterings_csv_dsn, index_col=0  ) .set_index ( ['sample_ID'] ) - 1
list_columns_str = df_cell_clusters_in.columns.values.tolist()
list_columns_int = [ int(c) for c in list_columns_str ]
rn_dict = dict ( zip ( list_columns_str, list_columns_int ) )
df_cell_clusters = df_cell_clusters_in.rename ( columns = rn_dict )
plog ( logfile, '\n\n df_cell_clusters:\n\n' , df_cell_clusters )



result =  MLE_cell_cluster_models ( df_UMI_counts, df_cell_clusters )
  
dict_MLE_results = result[0]
df_NLL = result[1] 
  
plog ( logfile, '\n\n df_NLL:\n\n' , df_NLL )  
  
  
  
  
f = open( dict_MLE_results_dsn, 'wb' )    
pickle.dump( dict_MLE_results, f)           
f.close()    

df_NLL.to_pickle ( NLL_LR_stat_summary_dsn )
df_cell_clusters.to_pickle ( Mclust_clusterings_pkl_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end  program  MLE_cell_cluster_models.py                                             # 
#                                                                                       #   
#########################################################################################
#########################################################################################



