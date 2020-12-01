

#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  MLE_cell_clusters_group_id.py                                         # 
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

logfile_txt =  "MLE_cell_clusters_group_id B cells.txt"

dict_MLE_results_pkl = "MLE_cell_clusters_group_id B cells.pkl"
# NLL_LR_stat_summary_pkl = "NLL_cell_clusters_group_id B cells.pkl"

column_data_select_pkl = "sce_column_data B cells.pkl"
counts_select_pkl = "sce_counts_B_cells.pkl"


data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle outputs
dict_MLE_results_dsn = data_path / dict_MLE_results_pkl
# NLL_LR_stat_summary_dsn = data_path / NLL_LR_stat_summary_pkl


# pickle inputs
column_data_select_dsn  = data_path / column_data_select_pkl 
counts_select_dsn = data_path / counts_select_pkl

######################################################################################################################################

df_column_data_select = pd.read_pickle ( column_data_select_dsn )
df_counts_select = pd.read_pickle ( counts_select_dsn )

plog ( logfile, '\n\n df_column_data_select:\n\n' , df_column_data_select )
plog ( logfile, '\n\n df_counts_select:\n\n' , df_counts_select )


n_clusters = 1 + df_column_data_select['cluster_group_id'].max()
df_cell_clusters = df_column_data_select[['cluster_group_id']].rename ( columns = {'cluster_group_id': n_clusters} )
plog ( logfile, '\n\n df_cell_clusters:\n\n' , df_cell_clusters )



result =  MLE_cell_cluster_models ( df_counts_select, df_cell_clusters )
  
dict_MLE_results = result[0]
df_NLL = result[1] 
  
plog ( logfile, '\n\n df_NLL:\n\n' , df_NLL )  
  
  
  
  
f = open( dict_MLE_results_dsn, 'wb' )    
pickle.dump( dict_MLE_results, f)           
f.close()    

# df_NLL.to_pickle ( NLL_LR_stat_summary_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  MLE_cell_clusters_group_id.py                                           # 
#                                                                                       #   
#########################################################################################
#########################################################################################



