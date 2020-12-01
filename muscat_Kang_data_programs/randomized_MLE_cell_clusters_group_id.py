

#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  randomized_MLE_cell_clusters_group_id.py                              # 
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

logfile_txt =  "randomized_MLE_cell_clusters_group_id B cells_11_randomizations.txt"

dict_MLE_results_list_pkl = "randomized_MLE_cell_clusters_group_id B cells_list_11_randomizations.pkl"
randomized_clusters_pkl = "randomized_cell_clusters_group_id B cells_11_randomizations.pkl"


column_data_select_pkl = "sce_column_data B cells.pkl"
counts_select_pkl = "sce_counts_B_cells.pkl"


data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle outputs
dict_MLE_results_list_dsn = data_path / dict_MLE_results_list_pkl
randomized_clusters_dsn = data_path / randomized_clusters_pkl
# NLL_LR_stat_summary_dsn = data_path / NLL_LR_stat_summary_pkl


# pickle inputs
column_data_select_dsn  = data_path / column_data_select_pkl 
counts_select_dsn = data_path / counts_select_pkl

######################################################################################################################################

n_randomizations = 11


df_column_data_select = pd.read_pickle ( column_data_select_dsn )
df_counts_select = pd.read_pickle ( counts_select_dsn )

plog ( logfile, '\n\n df_column_data_select:\n\n' , df_column_data_select )
plog ( logfile, '\n\n df_counts_select:\n\n' , df_counts_select )


n_clusters = 1 + df_column_data_select['cluster_group_id'].max()
df_cell_clusters = df_column_data_select[['cluster_group_id']]
plog ( logfile, '\n\n df_cell_clusters:\n\n' , df_cell_clusters )


cluster_list = df_cell_clusters['cluster_group_id'].values.tolist()
df_randomized_clusters_list = []

for randomization in range ( n_randomizations ):
  print ( 'randomization ', randomization )
  np.random.shuffle( cluster_list )
  df_randomized_cl = pd.DataFrame ( index=df_cell_clusters.index, data=cluster_list, columns=[ randomization ] )  
  df_randomized_clusters_list.append ( df_randomized_cl )
df_randomized_clusters = pd.concat ( df_randomized_clusters_list, axis=1, sort=False )
plog ( logfile, '\n\n df_randomized_clusters: \n\n', df_randomized_clusters )
pdline ( logfile )
  


dict_MLE_results_list = []
for randomization in range ( n_randomizations ):
  print ( 'randomization ', randomization )
  plog ( logfile, '\n\n\n randomization: ', randomization )
  df_cell_clusters_parm = df_randomized_clusters[[randomization]].rename ( columns = {randomization: n_clusters} )
  result =  MLE_cell_cluster_models ( df_counts_select, df_cell_clusters_parm )
  dict_MLE_results = result[0]
  dict_MLE_results_list.append ( dict_MLE_results )
  df_NLL = result[1]   
  plog ( logfile, '\n\n df_NLL:\n\n' , df_NLL )  
    
  
  
  
f = open( dict_MLE_results_list_dsn, 'wb' )    
pickle.dump( dict_MLE_results_list, f)           
f.close()    

df_randomized_clusters.to_pickle ( randomized_clusters_dsn )

# df_NLL.to_pickle ( NLL_LR_stat_summary_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  randomized_MLE_cell_clusters_group_id.py                                # 
#                                                                                       #   
#########################################################################################
#########################################################################################



