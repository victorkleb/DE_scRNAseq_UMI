

# https://stackoverflow.com/questions/5284646/rank-items-in-an-array-using-python-numpy-without-sorting-array-twice





#########################################################################################
#########################################################################################
#                                                                                       #    
#  start program  compare_z_scores_to_muscat_group_id.py                                # 
#                                                                                       #   
#########################################################################################

import pandas as pd
import numpy  as np

from scipy.stats import rankdata

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
 
 

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

def plot_consistency ( base_stat ):

  df_ranks_sort_base_stat = df_plot_stats_ranks.sort_values ( [ base_stat ] )

  base_rank_list =   df_ranks_sort_base_stat[base_stat].values.astype (np.int).tolist() 
  compare_list = stat_list.copy()
  compare_list.remove ( base_stat )

  df_consistent_fraction = df_ranks_sort_base_stat[[ base_stat ]]

  for stat in compare_list: 
    compare_ranks_array =   df_ranks_sort_base_stat[ stat ].values
    consistent_fraction_list = [ ( np.sum ( compare_ranks_array[:i] <=i ) ) / ( 1 + i ) for i in base_rank_list ]
    df_consistent_fraction [ stat ] = consistent_fraction_list


  title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells \n compare Differential Expression - use gene ranks to visualize CONSISTENCY between methods'
  title2 = '\n muscat pseudo bulk statistics and z-scores from  pi-hat range analysis' 
  fig, ax1 = plt.subplots()  
  titl = title1 + title2  
  plt.title( titl, fontsize=6.5 )
  
  for stat in compare_list:   
    ax1.plot (  df_consistent_fraction[base_stat],	df_consistent_fraction[stat],  c=color_dict[ stat ], linewidth=0.5,  label=stat )

  ax1.legend ( loc = 'lower right' , prop={'size': 5})
  ax1.set_xlabel ( base_stat + ' descending rank \n highest Differential Expression is on left', fontsize=6 )	
  ax1.set_ylabel ( 'fraction of genes with equal or higher \n Differential Expression rank' , fontsize=6 )	
  ax1.tick_params(labelsize=5 )           
  pdf_pages.savefig(fig, transparent=True )


  

def plot_inconsistency ( base_stat ):

  df_ranks_sort_base_stat = df_plot_stats_ranks.sort_values ( [ base_stat ] )

  base_rank_list =   df_ranks_sort_base_stat[base_stat].values.astype (np.int).tolist() 
  compare_list = stat_list.copy()
  compare_list.remove ( base_stat )

  df_inconsistent_ranks = df_ranks_sort_base_stat[[ base_stat ]]

  for stat in compare_list: 
    compare_ranks_array =   df_ranks_sort_base_stat[ stat ].values
    inconsistent_rank_list = [  np.max ( compare_ranks_array[:i] )  for i in base_rank_list ]
    df_inconsistent_ranks [ stat ] = inconsistent_rank_list

	
  title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells \n compare Differential Expression - use gene ranks to visualize INCONSISTENCY between methods'  
  title2 = '\n muscat pseudo bulk statistics and z-scores from  pi-hat range analysis' 
  fig, ax1 = plt.subplots()  
  titl = title1 + title2  
  plt.title( titl, fontsize=6.5 )
  
  for stat in compare_list:   
    ax1.plot (  df_inconsistent_ranks[base_stat],	df_inconsistent_ranks[stat],  c=color_dict[ stat ], linewidth=0.5,  label=stat )

  ax1.legend ( loc = 'lower right' , prop={'size': 5})
  ax1.set_xlabel ( base_stat + ' descending rank \n highest Differential Expression is on left', fontsize=6 )	
  ax1.set_ylabel ( 'maximum rank in alternate methods' , fontsize=6 )	
  ax1.tick_params(labelsize=5 )           
  pdf_pages.savefig(fig, transparent=True )


  
  

pctl_list =  [ .001, .005, .01, .05, .1, .25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ]

#########################################################################################  
 
logfile_txt = "compare_z_scores_to_muscat_group_id B cells.txt"
plot_pdf = "compare_z_scores_to_muscat_group_id B cells.pdf"

df_pi_hat_z_scores_41_pkl =  "pi_hat_z_scores_true_v_randomized_cell_clusters_group_id B cells_41_randomizations.pkl"
muscat_edgeR_csv = "edgeR_DE_B_cells.csv"
muscat_DESeq2_csv = "DESeq2_DE_B_cells.csv"
muscat_limma_trend_csv = "limma_trend_DE_B_cells.csv"
muscat_limma_voom_csv = "limma_voom_DE_B_cells.csv"


data_path = Path ( r"D:/scRNA-seq_DE/Crowell_Kang18_8vs8" )



####  log output
logfile_dsn  =  data_path / logfile_txt
logfile = open ( logfile_dsn,'w')

##  PDF output
plot_dsn = data_path /  plot_pdf
pdf_pages = PdfPages( plot_dsn )



# pickle input 
df_pi_hat_z_scores_dsn = data_path / df_pi_hat_z_scores_41_pkl

# csv input
muscat_edgeR_dsn = data_path / muscat_edgeR_csv
muscat_DESeq2_dsn =  data_path / muscat_DESeq2_csv
muscat_limma_trend_dsn = data_path / muscat_limma_trend_csv
muscat_limma_voom_dsn =  data_path / muscat_limma_voom_csv

######################################################################################################################################

color_dict = { 'z_score':'black', 'edgeR_F':'red',  'abs_DESeq2_stat':'blue', 'abs_limma_trend_t':'cyan', 'abs_limma_voom_t':'violet' }



df_z_scores = pd.read_pickle ( df_pi_hat_z_scores_dsn )
plog ( logfile, '\n\n\n\n df_z_scores: \n\n', df_z_scores )
plog ( logfile, '\n\n\n\n distribution of z_scores and p-values: \n', df_z_scores[['z_score', 'p_value', 'adj_p_value']] .describe ( percentiles = pctl_list ) )


df_muscat_edgeR = pd.read_csv ( muscat_edgeR_dsn )[['gene', 'F','p_val', 'p_adj.loc']].set_index ( ['gene'] ).rename ( columns={ 'p_val':'edgeR_p_val', 'p_adj.loc':'edgeR_p_adj', 'F':'edgeR_F'} )
plog ( logfile, '\n\n\n\n df_muscat_edgeR: \n\n', df_muscat_edgeR )

df_muscat_DESeq2 = pd.read_csv ( muscat_DESeq2_dsn )[['gene', 'stat','p_val', 'p_adj.loc']].set_index ( ['gene'] ).rename ( columns={ 'p_val':'DESeq2_p_val', 'p_adj.loc':'DESeq2_p_adj', 'stat':'DESeq2_stat'} )
df_muscat_DESeq2['abs_DESeq2_stat'] = np.abs ( df_muscat_DESeq2['DESeq2_stat'] )
plog ( logfile, '\n\n\n\n df_muscat_DESeq2: \n\n', df_muscat_DESeq2 )

df_muscat_limma_trend = pd.read_csv ( muscat_limma_trend_dsn )[['gene', 't','p_val', 'p_adj.loc']].set_index ( ['gene'] ).rename ( columns={ 'p_val':'limma_trend_p_val', 'p_adj.loc':'limma_trend_p_adj', 't':'limma_trend_t'} )
df_muscat_limma_trend['abs_limma_trend_t'] = np.abs ( df_muscat_limma_trend['limma_trend_t'] )
plog ( logfile, '\n\n\n\n df_muscat_limma_trend: \n\n', df_muscat_limma_trend )

df_muscat_limma_voom = pd.read_csv ( muscat_limma_voom_dsn )[['gene', 't','p_val', 'p_adj.loc']].set_index ( ['gene'] ).rename ( columns={ 'p_val':'limma_voom_p_val', 'p_adj.loc':'limma_voom_p_adj', 't':'limma_voom_t'} )
df_muscat_limma_voom['abs_limma_voom_t'] = np.abs ( df_muscat_limma_voom['limma_voom_t'] )
plog ( logfile, '\n\n\n\n df_muscat_limma_voom: \n\n', df_muscat_limma_voom )

df_muscat_stats = pd.concat ( [ df_muscat_edgeR, df_muscat_DESeq2, df_muscat_limma_trend, df_muscat_limma_voom ], axis=1, sort=False )
plog ( logfile, '\n\n\n\n df_muscat_stats: \n\n', df_muscat_stats )



df_z_scores['muscat_match'] = False
df_z_scores['muscat_match'].loc [ df_muscat_edgeR.index ] = True
plog ( logfile, '\n\n\n\n df_z_scores[muscat_match].sum(): \n\n', df_z_scores['muscat_match'].sum() )

plog ( logfile, '\n\n\n\n restrict:  z_scores match to muscat pseudo-bulk output \n distribution of z_scores and p-values: \n', \
df_z_scores[['z_score', 'p_value']].loc [  df_z_scores['muscat_match'] ].describe ( percentiles = pctl_list ) )
plog ( logfile, '\n\n\n\n restrict:  z_scores NO match to muscat pseudo-bulk output \n distribution of z_scores and p-values: \n', \
df_z_scores[['z_score', 'p_value']].loc [ ~ df_z_scores['muscat_match'] ].describe ( percentiles = pctl_list ) )

df_z_scores_muscat_match = df_z_scores.merge ( df_muscat_stats, how='inner', left_index=True, right_index=True )
df_z_scores_muscat_match['log10_edgeR_F'] = np.log10 ( df_z_scores_muscat_match['edgeR_F'] )
plog ( logfile, '\n\n\n\n df_z_scores_muscat_match: \n\n', df_z_scores_muscat_match )


df_corr = df_z_scores_muscat_match[['p_value', 'edgeR_p_val', 'DESeq2_p_val', 'limma_trend_p_val', 'limma_voom_p_val', 'z_score', 'edgeR_F',  'abs_DESeq2_stat', 'abs_limma_trend_t', 'abs_limma_voom_t']].corr( method='spearman' )
plog ( logfile, '\n\n\n\n df_corr: \n\n', df_corr )
     
df_corr = df_z_scores_muscat_match[['p_value', 'edgeR_p_val', 'DESeq2_p_val', 'limma_trend_p_val', 'limma_voom_p_val' ]].corr( method='spearman' )
plog ( logfile, '\n\n\n\n df_corr: \n\n', df_corr )
pdline ( logfile )


stat_list = ['z_score', 'edgeR_F',  'abs_DESeq2_stat', 'abs_limma_trend_t', 'abs_limma_voom_t']
df_plot_stats_ranks =  df_z_scores_muscat_match[ stat_list ].rank ( ascending=False )





title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare pi-hat range analysis to muscat pseudo bulk \n edgeR F statistic vs z-scores ' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['z_score'],	df_z_scores_muscat_match['log10_edgeR_F'],  c='k', s=1  )
ax1.set_xlabel ( 'z_score', fontsize=7 )	
ax1.set_ylabel ( 'log10 (edgeR_F)' , fontsize=7 )	
ax1.tick_params(labelsize=6)           
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare pi-hat range analysis to muscat pseudo bulk \n DESeq2_stat vs z-scores ' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['z_score'],	df_z_scores_muscat_match['abs_DESeq2_stat'],  c='k', s=1  )
ax1.set_xlabel ( 'z_score', fontsize=7 )	
ax1.set_ylabel ( 'abs_DESeq2_stat' , fontsize=7 )	
ax1.tick_params(labelsize=6)           
pdf_pages.savefig(fig, transparent=True )



title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare pi-hat range analysis to muscat pseudo bulk \n limma_trend_t vs z-scores ' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['z_score'],	df_z_scores_muscat_match['abs_limma_trend_t'],  c='k', s=1  )
ax1.set_xlabel ( 'z_score', fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_trend_t' , fontsize=7 )	
ax1.tick_params(labelsize=6)         
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare pi-hat range analysis to muscat pseudo bulk \n abs_limma_voom_t vs z-scores ' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['z_score'],	df_z_scores_muscat_match['abs_limma_voom_t'],  c='k', s=1  )
ax1.set_xlabel ( 'z_score', fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_voom_t' , fontsize=7 )	
ax1.tick_params(labelsize=6)           
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_DESeq2_stat vs edgeR F statistic' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['log10_edgeR_F'],	df_z_scores_muscat_match['abs_DESeq2_stat'],  c='k', s=1  )
ax1.set_xlabel ( 'log10 (edgeR_F)' , fontsize=7 )	
ax1.set_ylabel ( 'abs_DESeq2_stat', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_limma_trend_t vs edgeR F statistic' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['log10_edgeR_F'],	df_z_scores_muscat_match['abs_limma_trend_t'],  c='k', s=1  )
ax1.set_xlabel ( 'log10 (edgeR_F)' , fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_trend_t', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_limma_voom_t vs edgeR F statistic' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['log10_edgeR_F'],	df_z_scores_muscat_match['abs_limma_voom_t'],  c='k', s=1  )
ax1.set_xlabel ( 'log10 (edgeR_F)' , fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_voom_t', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_limma_trend_t vs abs_DESeq2_stat' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['abs_DESeq2_stat'],	df_z_scores_muscat_match['abs_limma_trend_t'],  c='k', s=1  )
ax1.set_xlabel ( 'abs_DESeq2_stat' , fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_trend_t', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_limma_voom_t vs abs_DESeq2_stat' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['abs_DESeq2_stat'],	df_z_scores_muscat_match['abs_limma_voom_t'],  c='k', s=1  )
ax1.set_xlabel ( 'abs_DESeq2_stat' , fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_voom_t', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


title1 = 'Crowell_Kang18_8vs8 scRNA-seq data - B cells: \n compare muscat pseudo bulk \n abs_limma_voom_t vs abs_limma_trend_t' 
fig, ax1 = plt.subplots()  
titl = title1  
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_z_scores_muscat_match['abs_limma_trend_t'],	df_z_scores_muscat_match['abs_limma_voom_t'],  c='k', s=1  )
ax1.set_xlabel ( 'abs_limma_trend_t' , fontsize=7 )	
ax1.set_ylabel ( 'abs_limma_voom_t', fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )




for base_stat in stat_list:  
  plot_consistency  ( base_stat )

  
for base_stat in stat_list: 
  plot_inconsistency ( base_stat )


pdf_pages.close()  
 


logfile.close()

#########################################################################################
#                                                                                       #    
#  end program  compare_z_scores_to_muscat_group_id.py                                  # 
#                                                                                       #   
#########################################################################################
#########################################################################################



