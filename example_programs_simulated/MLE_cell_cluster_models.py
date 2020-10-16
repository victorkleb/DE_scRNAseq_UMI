



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

logfile_txt =  "MLE_cell_cluster_models.txt"
dict_MLE_results_pkl = "MLE_cell_clusters_all_genes.pkl"
NLL_LR_stat_summary_pkl = "NLL_LR_stat_summary_all_genes.pkl"

synthetic_UMI_data_pkl = "synthetic_UMI_data.pkl"
clusters_pkl = "cell_clusters.pkl"


data_folder = Path ( r"D:/scRNA-seq_DE/simulated" )
data_path = Path ( data_folder )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle outputs
dict_MLE_results_dsn = data_path / dict_MLE_results_pkl
NLL_LR_stat_summary_dsn = data_path / NLL_LR_stat_summary_pkl


# pickle inputs
synthetic_UMI_data_dsn  =  data_path / synthetic_UMI_data_pkl
clusters_dsn = data_path / clusters_pkl 

######################################################################################################################################

df_UMI_counts = pd.read_pickle ( synthetic_UMI_data_dsn )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )

df_cell_clusters =  pd.read_pickle ( clusters_dsn )
plog ( logfile, '\n\n df_cell_clusters:\n\n' , df_cell_clusters )



result =  MLE_cell_cluster_models ( df_UMI_counts, df_cell_clusters )
  
dict_MLE_results = result[0]
df_NLL = result[1] 

plog ( logfile, '\n\n df_NLL:\n\n' , df_NLL )  
  
    
  
f = open( dict_MLE_results_dsn, 'wb' )    
pickle.dump( dict_MLE_results, f)           
f.close()    

df_NLL.to_pickle ( NLL_LR_stat_summary_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end  program  MLE_cell_cluster_models.py                                             # 
#                                                                                       #   
#########################################################################################
#########################################################################################



