



#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  prepare_synthetic_data.py                                             # 
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


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]
 
########################################################################################

logfile_txt =  "prepare_synthetic_data.txt"
synthetic_UMI_data_pkl = "synthetic_UMI_data.pkl"
clusters_pkl = "cell_clusters.pkl"

data_folder = r"D:/scRNA-seq_DE/simulated"
data_path = Path ( data_folder )


####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

####  pickle outputs
synthetic_UMI_data_dsn  =  data_path / synthetic_UMI_data_pkl
clusters_dsn = data_path / clusters_pkl 



######################################################################################################################################

n_cells_per_cluster = 500 
n_genes = 10000
background_parm  = 0.0002
DE_parm = 0.005
Poisson_rate_factor = 0.0001


df_UMI_counts_cluster_1 = pd.DataFrame ( index=range ( n_genes ) ) 
df_UMI_counts_cluster_2 = pd.DataFrame ( index=range ( n_genes ) ) 


arr_log_total_counts =   np.arange ( 4.5, 6.0,  1.5/n_cells_per_cluster ) 
arr_total_counts = np.floor ( 10** arr_log_total_counts )  
plog ( logfile, '\n\n distributon of arr_log_total_counts: \n\n', desc_list ( arr_total_counts ) )


background_exp_param = background_parm + np.zeros ( n_genes )  

DE_exp_parm_1 = np.copy ( background_exp_param )
DE_exp_parm_1[0:200] = DE_parm

DE_exp_parm_2 = np.copy ( background_exp_param )
DE_exp_parm_2[200:400] = DE_parm




plog  ( logfile,  '\n\n distribution for cells with Differentially Expressed genes 0-199 \n\n')

for i in range ( n_cells_per_cluster ):
  arr_Poisson_lambda = Poisson_rate_factor * DE_exp_parm_1 * arr_total_counts [i]
  arr_poisson_random = np.random.poisson ( arr_Poisson_lambda )
  df_UMI_counts_cluster_1[i] = arr_poisson_random  

plog ( logfile, '\n\n df_UMI_counts_cluster_1:\n\n' , df_UMI_counts_cluster_1 )

ser_gene_total_counts = df_UMI_counts_cluster_1.sum( axis=1 )
plog ( logfile, '\n\n ser_gene_total_counts:\n\n' , ser_gene_total_counts )

ser_cell_total_counts = df_UMI_counts_cluster_1.sum()
plog ( logfile, '\n\n ser_cell_total_counts:\n\n' , ser_cell_total_counts )

Log10_ser_gene_total_counts = np.log10 ( 0.001 + ser_gene_total_counts )
Log10_ser_cell_total_counts = np.log10 ( ser_cell_total_counts )


plog ( logfile, '\n\n distribution of total cell counts :\n\n' , ser_cell_total_counts.describe ( percentiles = pctl_list ) )
plog ( logfile, '\n\n distribution of total gene counts :\n\n' , ser_gene_total_counts.describe ( percentiles = pctl_list ) )

plog ( logfile, '\n\n distribution of Log10 ( total cell counts ) :\n\n' , Log10_ser_cell_total_counts.describe ( percentiles = pctl_list ) )  
plog ( logfile, '\n\n distribution of Log10 ( total gene counts ) :\n\n' , Log10_ser_gene_total_counts.describe ( percentiles = pctl_list ) )
 
pdline ( logfile )  
  
  
  
plog  ( logfile,  '\n\n distribution for cells with Differentially Expressed genes 200-399 \n\n')

for i in range ( n_cells_per_cluster ):
  arr_Poisson_lambda = Poisson_rate_factor * DE_exp_parm_2 * arr_total_counts [i]
  arr_poisson_random = np.random.poisson ( arr_Poisson_lambda )
  df_UMI_counts_cluster_2[ i + n_cells_per_cluster ] = arr_poisson_random  

plog ( logfile, '\n\n df_UMI_counts_cluster_2:\n\n' , df_UMI_counts_cluster_2 )

ser_gene_total_counts = df_UMI_counts_cluster_2.sum( axis=1 )
plog ( logfile, '\n\n ser_gene_total_counts:\n\n' , ser_gene_total_counts )

ser_cell_total_counts = df_UMI_counts_cluster_2.sum()
plog ( logfile, '\n\n ser_cell_total_counts:\n\n' , ser_cell_total_counts )

Log10_ser_gene_total_counts = np.log10 ( 0.001 + ser_gene_total_counts )
Log10_ser_cell_total_counts = np.log10 ( ser_cell_total_counts )
  
plog ( logfile, '\n\n distribution of total cell counts :\n\n' , ser_cell_total_counts.describe ( percentiles = pctl_list ) )
plog ( logfile, '\n\n distribution of total gene counts :\n\n' , ser_gene_total_counts.describe ( percentiles = pctl_list ) )

plog ( logfile, '\n\n distribution of Log10 ( total cell counts ) :\n\n' , Log10_ser_cell_total_counts.describe ( percentiles = pctl_list ) )  
plog ( logfile, '\n\n distribution of Log10 ( total gene counts ) :\n\n' , Log10_ser_gene_total_counts.describe ( percentiles = pctl_list ) )
 
pdline ( logfile ) 
  
  
  
df_UMI_counts = pd.concat ( [ df_UMI_counts_cluster_1, df_UMI_counts_cluster_2 ], axis=1, sort=False )    


plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )

ser_gene_total_counts = df_UMI_counts.sum( axis=1 )
plog ( logfile, '\n\n ser_gene_total_counts:\n\n' , ser_gene_total_counts )

ser_cell_total_counts = df_UMI_counts.sum()
plog ( logfile, '\n\n ser_cell_total_counts:\n\n' , ser_cell_total_counts )

Log10_ser_gene_total_counts = np.log10 ( 0.001 + ser_gene_total_counts )
Log10_ser_cell_total_counts = np.log10 ( ser_cell_total_counts )


plog ( logfile, '\n\n distribution of total cell counts :\n\n' , ser_cell_total_counts.describe ( percentiles = pctl_list ) )
plog ( logfile, '\n\n distribution of total gene counts :\n\n' , ser_gene_total_counts.describe ( percentiles = pctl_list ) )

plog ( logfile, '\n\n distribution of Log10 ( total cell counts ) :\n\n' , Log10_ser_cell_total_counts.describe ( percentiles = pctl_list ) )  
plog ( logfile, '\n\n distribution of Log10 ( total gene counts ) :\n\n' , Log10_ser_gene_total_counts.describe ( percentiles = pctl_list ) )
 

df_zeros = ( df_UMI_counts < 0.5 )
plog ( logfile, '\n\n df_zeros:\n\n' , df_zeros )
plog ( logfile, '\n\n df_zeros.sum.sum :\n\n' , df_zeros.sum().sum()  )
 
df_non_zeros = ( df_UMI_counts > 0.5 ) 
plog ( logfile, '\n\n df_non_zeros:\n\n' , df_non_zeros )
plog ( logfile, '\n\n df_non_zeros.sum.sum :\n\n' , df_non_zeros.sum().sum()  ) 


 
 

list_clusters = n_cells_per_cluster * [0] + n_cells_per_cluster * [1]

df_clusters = pd.DataFrame ( data=list_clusters, columns=[2] )
plog ( logfile,  '\n\n df_clusters: \n\n', df_clusters ) 



df_clusters.to_pickle ( clusters_dsn )

df_UMI_counts.to_pickle ( synthetic_UMI_data_dsn )


logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  prepare_synthetic_data.py                                               # 
#                                                                                       #   
#########################################################################################
#########################################################################################



