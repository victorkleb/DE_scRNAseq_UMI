################################################################################################################################ 
##                                                                                                                            ##
##  start program: compute_binomial_deviance.py                                                                               ##
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION compute_binomial_deviance__genuine_and_randomized_counts
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

pctl_list =  [ .1, .25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ]

#########################################################################################  

logfile_txt = "compute_binomial_deviance.txt"

genuine_data_pkl = 'binomial deviance - genuine data.pkl'
synthetic_data_pkl = 'binomial deviance - randomized data.pkl'
summary_pkl = 'binomial deviance summary table.pkl' 

UMI_counts_csv = "sce_full_Zhengmix8eq_counts.csv"


data_path = Path ( r"D:/scRNA-seq_DE/Zhengmix8eq" )
 
 
####  log output
logfile_dsn  =  data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle outputs - pandas data frames
genuine_out_dsn  = data_path / genuine_data_pkl 
synthetic_out_dsn  = data_path / synthetic_data_pkl 
summary_out_dsn = data_path / summary_pkl 



# input
UMI_counts_dsn = data_path / UMI_counts_csv

######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )
pdline( logfile )



result_list = compute_binomial_deviance__genuine_and_randomized_counts ( df_UMI_counts )

df_binomial_deviance__genuine_data = result_list[0]
df_binomial_deviance__randomized_data = result_list[1]
df_binomial_deviance_summary_table = result_list[2]

plog ( logfile, '\n\n\n df_binomial_deviance__genuine_data: \n', df_binomial_deviance__genuine_data )
plog ( logfile, '\n\n\n df_binomial_deviance__randomized_data: \n', df_binomial_deviance__randomized_data )
plog ( logfile, '\n\n\n df_binomial_deviance_summary_table: \n', df_binomial_deviance_summary_table )


plog ( logfile, '\n\n\n distribution of df_binomial_deviance__genuine_data: \n', df_binomial_deviance__genuine_data.describe ( percentiles = pctl_list ) )
plog ( logfile, '\n\n\n distribution of df_binomial_deviance__randomized_data: \n', df_binomial_deviance__randomized_data.describe ( percentiles = pctl_list ) )



df_binomial_deviance__genuine_data.to_pickle ( genuine_out_dsn )
df_binomial_deviance__randomized_data.to_pickle ( synthetic_out_dsn )
df_binomial_deviance_summary_table.to_pickle ( summary_out_dsn )


logfile.close()
 
################################################################################################################################ 
##                                                                                                                            ##
##  end program: compute_binomial_deviance.py                                                                                 ##
##                                                                                                                            ## 
################################################################################################################################ 