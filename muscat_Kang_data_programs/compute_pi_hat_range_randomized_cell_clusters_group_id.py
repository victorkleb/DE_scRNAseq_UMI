

#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  compute_pi_hat_range_randomized_cell_clusters_group_id.py             # 
#                                                                                       #   
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



logfile_txt = "compute_pi_hat_range_randomized_cell_clusters_group_id B cells_11_randomizations.txt"
df_pi_hat_range_pkl =  "pi_hat_range_randomized_cell_clusters_group_id B cells_11_randomizations.pkl"


dict_MLE_results_list_pkl = "randomized_MLE_cell_clusters_group_id B cells_list_11_randomizations.pkl"



data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle output
df_pi_hat_range_dsn = data_path / df_pi_hat_range_pkl


# pickle input
dict_MLE_results_list_dsn = data_path / dict_MLE_results_list_pkl


######################################################################################################################################

group_id_clustering = 2

f = open( dict_MLE_results_list_dsn, 'rb' )    
dict_MLE_results_list = pickle.load( f)           
f.close()    

n_randomizations = len ( dict_MLE_results_list )



df_pi_hat_range_list = []

for randomization in  range ( n_randomizations ):
  print ( 'randomization ', randomization )
  plog ( logfile, '\n\n\n randomization: ', randomization )
  dict_MLE_results = dict_MLE_results_list [ randomization ]
  dict_cluster_results = dict_MLE_results[ group_id_clustering ] 
  df_pi_hat = dict_cluster_results[ 'pi_hat' ]
  
  list_pi_hat_columns = df_pi_hat.columns.values.tolist()  
  df_pi_hat['pi_hat_max'] = df_pi_hat[ list_pi_hat_columns ].max(axis=1)
  df_pi_hat['pi_hat_min'] = df_pi_hat[ list_pi_hat_columns ].min(axis=1)
  df_pi_hat['pi_hat_range'] = df_pi_hat['pi_hat_max'] - df_pi_hat['pi_hat_min']   
  df_pi_hat_range_list.append (  df_pi_hat[['pi_hat_range']].rename ( columns= {'pi_hat_range': randomization} ) )
  

df_pi_hat_range = pd.concat ( df_pi_hat_range_list, axis=1, sort=False )  
plog ( logfile, '\n\n\n df_pi_hat_range: \n\n', df_pi_hat_range )
plog ( logfile, '\n\n distribution of pi_hat_range: \n', df_pi_hat_range.describe ( percentiles = pctl_list ) )


df_pi_hat_range.to_pickle ( df_pi_hat_range_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  compute_pi_hat_range_randomized_cell_clusters_group_id.py               # 
#                                                                                       #   
#########################################################################################
#########################################################################################



