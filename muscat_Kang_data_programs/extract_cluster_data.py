################################################################################################################################ 
##                                                                                                                            ##
##  start program: extract_cluster_data.py                                                                                    ##
##                                                                                                                            ## 
################################################################################################################################


import pandas as pd
import numpy  as np



from scipy.stats import chi2_contingency
 
from pathlib import Path
import os


import sys
sys.path.append("D:/scRNA-seq_DE/DE_scRNAseq")


from  utilities import *
from FUNCTIONS_DE_scRNAseq  import *




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

def pv_table (  dsin0, col1,   col2 ):
  # print( 'col1: ', col1, '    col2: ', col2 )
  dsin = dsin0.copy()
  dsin ['count'] = 1
  pt = pd.pivot_table( dsin, values='count',  index=[ col1 ], columns=[ col2 ], aggfunc=np.sum, margins=True )
  ptfna = pt.fillna(0)
  pti = ptfna.astype(int)
  print ( '\n\n\n ',pti ,  file = logfile )
  g, p, dof, expctd = chi2_contingency( pti, lambda_="log-likelihood")
  print ( '\ng-statistic: ', format(g, '.2f') , '     p-value: ', format(p, '.6f'), '\n', file = logfile )
  
  return pti




pctl_list =  [ .1, .25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ]

#########################################################################################  

logfile_txt = "extract_cluster_data B cells.txt"

column_data_select_pkl = "sce_column_data B cells.pkl"
counts_select_pkl = "sce_counts_B_cells.pkl"

counts_csv = "sce_counts.csv"
column_data_csv = "sce_column_data.csv"


data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs 
column_data_select_dsn  = data_path / column_data_select_pkl 
counts_select_dsn = data_path / counts_select_pkl


# input
counts_dsn = data_path / counts_csv
column_data_dsn  = data_path / column_data_csv


######################################################################################################################################

df_column_data = pd.read_csv ( column_data_dsn, index_col=0 ).drop ( columns=['cell','sizeFactor', 'factor_ind'] ) .rename ( columns={'stim':'group_id'} )
plog ( logfile, '\n\n df_column_data:\n\n' , df_column_data )
plog ( logfile, '\n\n df_column_data - columns: \n\n', df_column_data.columns.values.tolist() ) 
plog ( logfile, '\n\n unique values of cluster_id: \n',  df_column_data['cluster_id'].unique().tolist() ) 


df_counts = pd.read_csv ( counts_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_counts:\n\n' , df_counts )
pdline( logfile )



plog ( logfile, '\n\n cross-tab numeric cluster vs alpha cluster_id: ' )
xtab = pv_table ( df_column_data,  'cluster', 'cluster_id' )


cluster_select =  'B cells'
df_column_data_select = df_column_data.loc [ df_column_data['cluster_id'] == cluster_select ]
plog ( logfile, '\n\n df_column_data_select:\n\n' , df_column_data_select )

list_group_id_values = df_column_data_select['group_id'].unique().tolist()
list_group_id_values.sort()

list_sample_id_values = df_column_data_select['sample_id'].unique().tolist()
list_sample_id_values.sort()

dict_rn_group_id = dict ( zip (  list_group_id_values, list ( range ( 1 + len ( list_group_id_values ) ) ) ) ) 
dict_rn_sample_id = dict ( zip (  list_sample_id_values, list ( range ( 1 + len ( list_sample_id_values ) ) ) ) ) 

df_column_data_select['cluster_group_id'] = df_column_data_select['group_id'].map ( dict_rn_group_id )
df_column_data_select['cluster_sample_id'] = df_column_data_select['sample_id'].map ( dict_rn_sample_id )

plog ( logfile, '\n\n df_column_data_select:\n\n' , df_column_data_select )



df_counts_select0 = df_counts [ df_column_data_select.index.values.tolist() ]
plog ( logfile, '\n\n df_counts_select0:\n\n' , df_counts_select0 )

ser_counts_select0_number_GT_0 = ( df_counts_select0 > 0 ).sum ( axis=1 )
plog ( logfile, '\n\n distribution of ser_counts_select0_number_GT_0:\n\n' , ser_counts_select0_number_GT_0.describe ( percentiles = pctl_list ) )

df_counts_select = df_counts_select0.loc [ ser_counts_select0_number_GT_0 >= 10 ]
plog ( logfile, '\n\n df_counts_select:\n\n' , df_counts_select )

df_column_data_select.to_pickle ( column_data_select_dsn )
df_counts_select.to_pickle ( counts_select_dsn )


logfile.close()
 
################################################################################################################################ 
##                                                                                                                            ##
##  ends program: extract_cluster_data.py                                                                                     ##
##                                                                                                                            ## 
################################################################################################################################ 