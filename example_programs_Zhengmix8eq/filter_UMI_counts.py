################################################################################################################################ 
################################################################################################################################ 
##                                                                                                                            ##
##  start program: filter_UMI_counts.py                                                                                       ##
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION filter_UMI_counts
#
################################################################################################################################

   
	

import pandas as pd
import numpy  as np

 
from pathlib import Path
import os


import sys
sys.path.append("D:/scRNA-seq_DE/DE_scRNAseq")


from  utilities import *
from FUNCTIONS_DE_scRNAseq  import *




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################
#########################################################################################  

n_genes = 195

logfile_txt = "filter_UMI_counts.txt"
filter_UMI_counts_pkl = 'filtered_UMI_counts.pkl'
filter_UMI_counts_csv = 'filtered_UMI_counts.csv'


binomial_deviance_pkl = 'binomial deviance - genuine data.pkl'
UMI_counts_csv = "sce_full_Zhengmix8eq_counts.csv"


data_path = Path ( r"D:/scRNA-seq_DE/Zhengmix8eq" )
 
 
####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')


#### pickle output  
filter_UMI_counts_pkl_dsn  = data_path / filter_UMI_counts_pkl 

# CSV output
filter_UMI_counts_csv_dsn  = data_path / filter_UMI_counts_csv



# CSV input
UMI_counts_dsn = data_path /  UMI_counts_csv

# pickle input
binomial_deviance_dsn = data_path / binomial_deviance_pkl

######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )
pdline( logfile )

df_binomial_deviance = pd.read_pickle ( binomial_deviance_dsn ) 


result = filter_UMI_counts ( df_UMI_counts, df_binomial_deviance, n_genes )
df_filtered_UMI_counts = result [ 0 ]
plog ( logfile, '\n\n df_filtered_UMI_counts:\n\n' , df_filtered_UMI_counts )


df_filtered_UMI_counts.to_pickle ( filter_UMI_counts_pkl_dsn )
df_filtered_UMI_counts.to_csv ( filter_UMI_counts_csv_dsn )


logfile.close()
 	

	
################################################################################################################################ 
##                                                                                                                            ##
##  end program: filter_UMI_counts.py                                                                                         ##
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################
