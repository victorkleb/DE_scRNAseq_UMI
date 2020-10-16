


#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  compare_gene_range_of_MLE_to_binomial_deviance.py                     # 
#                                                                                       #   
#########################################################################################


import pandas as pd
import numpy  as np


import matplotlib.pyplot as plt



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

logfile_txt =  "compare_gene_range_of_MLE_to_binomial_deviance.txt"
plot_png =  "compare_gene_range_of_MLE_to_binomial_deviance.png"

dict_MLE_results_pkl = "MLE_cell_clusters_all_genes.pkl"
binomial_deviance_saturated_all_genes_pkl = 'binomial deviance - genuine data.pkl'


data_folder = Path ( r"D:/scRNA-seq_DE/simulated" )
data_path = Path ( data_folder )


####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

###  PDF output
plot_dsn = data_path / plot_png

 

#### pickle inputs
dict_MLE_results_dsn = data_path / dict_MLE_results_pkl
binomial_deviance_saturated_all_genes_dsn  = data_folder / binomial_deviance_saturated_all_genes_pkl 

###########################################################################

f = open( dict_MLE_results_dsn, 'rb' )    
dict_MLE_results = pickle.load ( f )           
f.close()  


df_saturated_binomial_deviance =  pd.read_pickle ( binomial_deviance_saturated_all_genes_dsn ).rename ( columns={ 'genuine': 'saturated'} )
plog ( logfile, '\n\n df_saturated_binomial_deviance:\n' , df_saturated_binomial_deviance )

##################################

clustering = 2
plog ( logfile, '\n\n\n number of clusters: ', clustering ) 

  
results_dict = dict_MLE_results[ clustering ]
df_pi_hat = results_dict['pi_hat']  
plog ( logfile, '\n\n df_pi_hat:\n' , df_pi_hat )

df_pi_hat_sum = df_pi_hat.sum()
plog ( logfile, '\n\n df_pi_hat_sum.describe:\n' , df_pi_hat_sum.describe () )
  
  
  
df_analy = df_saturated_binomial_deviance.copy()
df_analy ['Log10_saturated_bin_dev'] = np.log10 ( df_analy['saturated'] )

df_analy['pi_hat_max'] = df_pi_hat.max(axis=1)
df_analy['pi_hat_min'] = df_pi_hat.min(axis=1)
df_analy['pi_hat_range'] = df_analy['pi_hat_max'] - df_analy['pi_hat_min'] 
df_analy['Log10_pi_hat_range'] =  np.log10 ( df_analy['pi_hat_range'] )

plog ( logfile, '\n\n df_analy:\n' , df_analy )
plog ( logfile, '\n\n df_analy.describe:\n' , df_analy[[ 'pi_hat_max', 'pi_hat_min', 'pi_hat_range', 'Log10_pi_hat_range', 'Log10_saturated_bin_dev' ]].describe ( percentiles = pctl_list ) )
  
df_analy_DE_genes = df_analy.loc [ df_analy.index < 400 ]

plog ( logfile, '\n\n df_analy_DE_genes:\n' , df_analy_DE_genes )
plog ( logfile, '\n\n df_analy_DE_genes.describe:\n' , df_analy_DE_genes[[ 'pi_hat_max', 'pi_hat_min', 'pi_hat_range', 'Log10_pi_hat_range','Log10_saturated_bin_dev' ]].describe ( percentiles = pctl_list ) )
  

  
title1 = 'SIMULATED data: 10,000 genes, 1,000 cells, 2 clusters'
title2 = "\n compare range of each gene's cluster MLEs to binomial deviance"
title3 = '\n Differentially Expressed genes are highlighted in red'
titl = title1 + title2 + title3

fig,ax = plt.subplots( figsize=(8.,9) )  

ax.set_title( titl, fontsize=12.0 )
  
ax.scatter ( df_analy ['Log10_saturated_bin_dev'], df_analy['Log10_pi_hat_range']  , s=0.1, color='black'  )
ax.scatter ( df_analy_DE_genes ['Log10_saturated_bin_dev'], df_analy_DE_genes['Log10_pi_hat_range']  , s=1.0, color='red'  )

ax.set_xlabel ( 'Log10 ( binomial deviance )', fontsize=10.5 )	
ax.set_ylabel ( 'Log10 ( MLE range )' , fontsize=10.5 )	
plt.xticks( fontsize= 9.5 )
plt.yticks( fontsize= 9.5 )


fig.savefig( plot_dsn, dpi=300 )  


  

logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  compare_gene_range_of_MLE_to_binomial_deviance.py                       # 
#                                                                                       #   
#########################################################################################
#########################################################################################



