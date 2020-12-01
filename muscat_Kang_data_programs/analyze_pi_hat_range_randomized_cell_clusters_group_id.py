
# https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.fdrcorrection.html
# https://stackoverflow.com/questions/25185205/calculating-adjusted-p-values-in-python


#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  analyze_pi_hat_range_randomized_cell_clusters_group_id.py             # 
#                                                                                       #   
#########################################################################################

import pandas as pd
import numpy  as np

import scipy.stats as st
 
import statsmodels.stats.multitest as multi

 
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

pctl_list =  [ .001, .005, .01, .05, .1, .25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ]

#########################################################################################  

logfile_txt = "analyze_pi_hat_range_randomized_cell_clusters_group_id B cells_41_randomizations.txt"

df_pi_hat_z_scores_pkl =  "pi_hat_z_scores_true_v_randomized_cell_clusters_group_id B cells_41_randomizations.pkl"

df_pi_hat_range_pkl =  "pi_hat_range_randomized_cell_clusters_group_id B cells_41_randomizations.pkl"
dict_MLE_results_pkl = "MLE_cell_clusters_group_id B cells.pkl"





data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle output 
df_pi_hat_z_scores_dsn = data_path / df_pi_hat_z_scores_pkl

# pickle inputs
df_pi_hat_range_dsn = data_path / df_pi_hat_range_pkl
dict_MLE_results_dsn = data_path / dict_MLE_results_pkl


######################################################################################################################################

group_id_clustering = 2


df_pi_hat_randomized_range = pd.read_pickle ( df_pi_hat_range_dsn )
plog ( logfile, '\n\n\n df_pi_hat_randomized_range: \n\n', df_pi_hat_randomized_range )  

df_corr_sp = df_pi_hat_randomized_range.corr ( method='spearman' )
plog ( logfile, '\n\n\n df_corr_sp: \n\n', df_corr_sp )  



f = open( dict_MLE_results_dsn, 'rb' )    
dict_MLE_results = pickle.load( f)           
f.close()    


dict_cluster_results = dict_MLE_results[ group_id_clustering ] 
df_pi_hat = dict_cluster_results[ 'pi_hat' ]
list_pi_hat_columns = df_pi_hat.columns.values.tolist()
df_pi_hat['max'] = df_pi_hat[ list_pi_hat_columns ].max(axis=1)
df_pi_hat['min'] = df_pi_hat[ list_pi_hat_columns ].min(axis=1)
df_pi_hat['range'] = df_pi_hat['max'] - df_pi_hat['min'] 
plog ( logfile, '\n\n\n df_pi_hat: \n\n', df_pi_hat )
plog ( logfile, '\n\n distribution of range: \n', df_pi_hat['range'].describe ( percentiles = pctl_list ) )


df_z_scores = df_pi_hat_randomized_range.mean ( axis= 1).to_frame ( name = 'mean' )
df_z_scores['std'] = df_pi_hat_randomized_range.std ( axis= 1)
df_z_scores['range'] = df_pi_hat['range']
df_z_scores['z_score'] = ( df_z_scores['range']  - df_z_scores['mean'] ) / df_z_scores['std']
df_z_scores['p_value'] = st.norm.sf(  df_z_scores['z_score'] ) 

arr_p_value = df_z_scores['p_value'].values
( arr_bool, arr_fdr ) =  multi.fdrcorrection( arr_p_value )
df_z_scores['adj_p_value'] = arr_fdr




plog ( logfile, '\n\n\n df_z_scores: \n\n', df_z_scores )
plog ( logfile, '\n\n distribution of z_scores and p-values: \n', df_z_scores[['z_score', 'p_value', 'adj_p_value']] .describe ( percentiles = pctl_list ) )


df_z_scores.to_pickle ( df_pi_hat_z_scores_dsn )



logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  analyze_pi_hat_range_randomized_cell_clusters_group_id.py               # 
#                                                                                       #   
#########################################################################################
#########################################################################################



