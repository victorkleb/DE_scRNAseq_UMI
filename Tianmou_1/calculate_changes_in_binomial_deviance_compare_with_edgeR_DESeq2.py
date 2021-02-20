


#########################################################################################
#                                                                                       #    
#  program  calculate_changes_in_binomial_deviance_compare_with_edgeR_DESeq2.py         # 
#                                                                                       #   
#########################################################################################
#                                                                                       #                    
#  invoke functions   MLE_cell_cluster_model                                            #
#                     change_in_binomial_deviance                                       #
#                                                                                       #   
#########################################################################################

import pandas as pd
import numpy  as np


from scipy.stats import f
from statsmodels.stats.multitest import  fdrcorrection 
  

from sklearn import metrics

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
 

  
 
from pathlib import Path
import os



import sys
sys.path.append("D:/scRNA-seq_DE/DE_scRNAseq")



pd.options.display.width = 180
pd.set_option('display.max_columns', 30)

#########################################################################################

def  MLE_cell_cluster_model  ( df_UMI_counts, ser_cell_clusters ):

  n_cells = df_UMI_counts.shape[1]
  n_genes = df_UMI_counts.shape[0]
  n_clusters = 1+ ser_cell_clusters.max()

  df_cell_cluster_gene_totals_list = []
  
  for cluster in range ( n_clusters ):
    cluster_cell_list = ser_cell_clusters.loc [ ser_cell_clusters==cluster ].index.values.tolist()
    df_UMI_counts_cluster = df_UMI_counts [ cluster_cell_list ]
    df_cell_cluster_gene_totals_list.append ( df_UMI_counts_cluster.sum ( axis=1 ).to_frame ( name = cluster ) )  
    
  df_cell_cluster_gene_total_counts = pd.concat ( df_cell_cluster_gene_totals_list, axis=1, sort=False ) 

  ser_cluster_total_counts = df_cell_cluster_gene_total_counts.sum() 
  df_pi_hat = df_cell_cluster_gene_total_counts. divide ( ser_cluster_total_counts )
  
  df_product_count_log_pi_hat = df_cell_cluster_gene_total_counts.multiply ( np.log ( 1.e-20 + df_pi_hat ) )
  NLL = - df_product_count_log_pi_hat.sum().sum()  
   
  dict_cluster_results =  { 'cell_cluster_gene_totals': df_cell_cluster_gene_total_counts, 'pi_hat':df_pi_hat, 'NLL': NLL  }
 
  return  dict_cluster_results 
 
 
  
def change_in_binomial_deviance ( df_base_pi_hat, df_nested_counts,  df_cell_data_select_base_nested ):
  
  arr_base_nested_segmentation = df_cell_data_select_base_nested.columns.values
  base_segmentation = arr_base_nested_segmentation[0]
  nested_segmentation = arr_base_nested_segmentation[1]  
  
  base_segmentation_list = df_cell_data_select_base_nested[ base_segmentation ].unique().tolist() 
  base_segmentation_list.sort()

  nested_ID_in_base_ID_list = []
  for cg in base_segmentation_list:
    nested_ID_list = df_cell_data_select_base_nested[ nested_segmentation ].loc [ df_cell_data_select_base_nested[ base_segmentation ] == cg ].unique().tolist()
    nested_ID_list.sort()
    nested_ID_in_base_ID_list.append ( nested_ID_list )  	
	
	
  
  df_rhs_list = []
  
  for base_segmentation in base_segmentation_list:
    arr_gen_counts = df_nested_counts[ nested_ID_in_base_ID_list[ base_segmentation] ].values 
    mat_gen_counts = np.asmatrix ( arr_gen_counts )
  
    pi_hat_vector = np.transpose ( np.asmatrix ( df_base_pi_hat[ base_segmentation ].values ) )  
    n_vector = mat_gen_counts.sum(axis=0)  
    n_pi_hat = np.array ( pi_hat_vector * n_vector )
    n_1_m_pi_hat =  np.array( n_vector ) - n_pi_hat 

    term1 = arr_gen_counts * np.log ( 1e-20 + np.nan_to_num( arr_gen_counts / ( 1e-20 + n_pi_hat ) ) )
    arr_ni_m_gen_counts = np.array( n_vector ) - arr_gen_counts
    term2 = arr_ni_m_gen_counts * np.log ( 1e-20 +  arr_ni_m_gen_counts / ( 1e-20 + n_1_m_pi_hat ) )
    rhs = 2 * ( term1 + term2 )  
  	
  
    df_rhs_list.append ( pd.DataFrame ( index=df_nested_counts.index, data = rhs ) )

  ser_change = pd.concat ( df_rhs_list, axis=1, sort=False ).sum ( axis=1 ) 
  return ser_change   
    
	

def calculate_changes_in_binomial_deviance ( df_counts ):
    
  dict_cluster_model = MLE_cell_cluster_model ( df_counts, df_cell_data['Group']  )
  print (  '\n\n\n cluster model - cell cluster gene totals: \n\n', dict_cluster_model[ 'cell_cluster_gene_totals' ] )
  cluster_tuple = ( 'Cluster model', dict_cluster_model['NLL'] )  
  
  dict_null_model = MLE_cell_cluster_model ( df_counts, df_cell_data['Null_model']  )
  print (  '\n\n\n null model - cell cluster gene totals: \n\n', dict_null_model[ 'cell_cluster_gene_totals' ] )
  null_tuple = ( 'Null model', dict_null_model['NLL'] )
      
  dict_saturated_model = MLE_cell_cluster_model ( df_counts, df_cell_data['Saturated_model']  )
  print (  '\n\n\n saturated model - cell cluster gene totals: \n\n', dict_saturated_model[ 'cell_cluster_gene_totals' ] )
  saturated_tuple = ( 'Saturated_model' , dict_saturated_model['NLL'] )
  
  NLL_tuple_list = [ null_tuple, cluster_tuple, saturated_tuple ]  
  df_NLL = pd.DataFrame ( data = NLL_tuple_list, columns=['model','NLL'] )
  print (  '\n\n\n df_NLL: \n\n', df_NLL )
 
 
  df_change_cluster_to_null = change_in_binomial_deviance ( dict_null_model['pi_hat'], dict_cluster_model['cell_cluster_gene_totals'],   df_cell_data[[ 'Null_model', 'Group' ]] ).to_frame ( name = 'cluster_to_null' ) 
  df_change_saturated_to_cluster = change_in_binomial_deviance ( dict_cluster_model['pi_hat'], dict_saturated_model['cell_cluster_gene_totals'],  df_cell_data[[ 'Group', 'Saturated_model' ]] ).to_frame ( name = 'saturated_to_cluster' )        
  df_changes = pd.concat ( [ df_change_cluster_to_null, df_change_saturated_to_cluster ],  axis=1, sort=False )
  df_changes['bin_dev_F'] = df_changes['cluster_to_null'] / df_changes['saturated_to_cluster'] * ( n_cells - 2 )
  df_changes['log10_bin_dev_F'] = np.log10 ( df_changes['bin_dev_F'] )
  df_changes['bin_dev_p_value'] =  df_changes['bin_dev_F'].apply ( lambda F: 1 - f.cdf( F, 1,  ( n_cells-1) ) )
  reject, padj_array = fdrcorrection( df_changes['bin_dev_p_value'].values, alpha=0.05  )
  df_changes['bin_dev_p_value_adj'] = padj_array  
  print (  '\n\n\n df_changes: \n\n', df_changes  ) 

  return ( df_changes, df_NLL )


  

def scatter_plot ( df_plot, x, y, x_label, y_label ):

  fig, ax1 = plt.subplots()  
  titl = title1 +  title2 + title3
  plt.title( titl, fontsize=5.5 )

  df_DE = df_plot.loc [ df_plot['DE']==1 ]
  df_not_DE = df_plot.loc [ df_plot['DE']==0 ]
 
  ax1.scatter ( df_DE[x], df_DE[y],  c='red', s=1, label='Differentially Expressed') 
  ax1.scatter ( df_not_DE[x], df_not_DE[y],  c='dodgerblue', s=0.1, label='NOT Differentially Expressed')  
  
  ax1.set_xlabel ( x_label, fontsize=6.1 )	
  ax1.set_ylabel ( y_label, fontsize=6.1 )	
  ax1.tick_params(labelsize=5.2) 

  ax1.legend ( loc='lower left', fontsize=5 )	
  pdf_pages.savefig(fig, transparent=True )
  
  


def ROC_plots ( ):
  global title1, title2, title3
  
  number_DE = df_compare_methods_notNA_DE['DE'].sum()  
  
  dict_method_ROCs = {}
  for stat,color in zip ( method_stat_list, line_color_list ):
    fpr, tpr, thresholds = metrics.roc_curve(  df_compare_methods_notNA_DE['DE'].values,  np.abs ( df_compare_methods_notNA_DE[stat].values ) ) 
    df_ROC = pd.DataFrame ( data = zip ( fpr,tpr ), columns = ['fpr','tpr'] )    
    AUC = metrics.roc_auc_score(  df_compare_methods_notNA_DE['DE'].values,  np.abs ( df_compare_methods_notNA_DE[stat].values ) ) 
    dict_method_ROCs[ stat ] = { 'df_ROC':df_ROC, 'AUC': AUC, 'color':color }       
 
    
  fig, ax1 = plt.subplots()
  titl = title1 + title2 + title3 
  plt.title( titl, fontsize=5.0 )

  AUC_list = [] 
  for stat in  method_stat_list:
    dict_ROC_AUC = dict_method_ROCs[stat]
    df_ROC = dict_ROC_AUC['df_ROC']
    AUC = dict_ROC_AUC['AUC']  
    AUC_list.append( AUC )	
    color = dict_ROC_AUC['color']  
   
    ax1.plot ( df_ROC['fpr'],	df_ROC['tpr'],   c=color, linewidth=1., label=stat + ' / ' + '{0:.3f}'.format(AUC)  )
    
  ax1.set_xlabel ( 'False positive rate', fontsize=6.1 )	
  ax1.set_ylabel ( 'True positive rate' , fontsize=6.1 )	
  ax1.tick_params(labelsize=5.2) 
  
  ax1.legend ( loc = 'lower right' , prop={'size': 5}, title='method/AUC', title_fontsize=5 )   
  pdf_pages.savefig(fig, transparent=True )  

  df_AUC = pd.DataFrame ( index = method_stat_list, data= AUC_list, columns=['AUC'] ) 
  return  df_AUC   
   
  
  
  
def plots ( expression_segment ):
  global df_compare_methods_notNA_DE, title1, title2, title3

  df_compare['log10_edgeR_F'] = np.log10 ( df_compare['edgeR_F'] )
  print (  '\n\n\n df_compare: \n\n', df_compare )  
   
  df_compare_methods_notNA = df_compare.loc [ df_compare['bin_dev_F'].notna() ]
  df_compare_methods_notNA_DE = df_compare_methods_notNA.merge ( df_trueDE, how='left', left_index=True, right_index=True )
  df_compare_methods_notNA_DE['DE'] = df_compare_methods_notNA_DE['DE'].fillna( 0 )    
  
  n_genes = df_compare.shape[0]  
  number_DE = df_compare_methods_notNA_DE['DE'].sum()
  
  
  title1 =  "Mou et al. Beta-Poisson simulated data with log fold-change level 1 - " + expression_segment + " Expressed genes"
 
  title2 = '\n compare DESeq2 statistic to edgeR F-statistic \n ' + str ( n_genes ) + ' genes'
  title3 = '\n (for DESeq2, added pseudo count of 1)'
  scatter_plot ( df_compare_methods_notNA_DE,  'log10_edgeR_F', 'DESeq2_stat', 'Log10 ( edgeR F-statistic )', 'DESeq2 statistic' )

  title2 = '\n compare edgeR F-statistic to binomial deviance F statistic \n ' + str ( n_genes ) + ' genes'
  title3 = ''
  scatter_plot ( df_compare_methods_notNA_DE,  'log10_bin_dev_F', 'log10_edgeR_F', 'Log10 ( binomial deviance F statistic )', 'Log10 ( edgeR F-statistic )' )

  title2 = '\n compare DESeq2 statistic to binomial deviance F statistic \n ' + str ( n_genes ) + ' genes'
  title3 = '\n (for DESeq2, added pseudo count of 1)' 
  scatter_plot ( df_compare_methods_notNA_DE,  'log10_bin_dev_F', 'DESeq2_stat', 'Log10 ( binomial deviance F statistic )', 'DESeq2 statistic' ) 
 
  title2 = '\n ROC curves for edgeR, DESeq2, and proposed method'
  title3 = '\n total genes: ' + str ( n_genes ) + '   number DE: ' +    '{0:.0f}'.format( number_DE )    
  df_AUC = ROC_plots()   
  
  return   df_AUC 
 
#########################################################################################
 
plot_pdf    = "compare_methods.pdf"
Excel_xlsx  = "compare_methods.xlsx"


trueDE_csv    =  "trueDE.csv"

counts_low_csv  = "counts_low.csv"
counts_high_csv = "counts_high.csv"

edgeR_DESeq2_low_csv  = "edgeR_and_DESeq_stats_low.csv"
edgeR_DESeq2_high_csv = "edgeR_and_DESeq_stats_high.csv"



work_folder =  "D:/scRNA-seq_DE/simulation_Tianmou_1"

data_path = Path ( work_folder )



##  PDF output
plot_dsn = data_path /  plot_pdf
pdf_pages = PdfPages( plot_dsn )

# Excel output
Excel_dsn = data_path / Excel_xlsx 
writer = pd.ExcelWriter( Excel_dsn, engine='xlsxwriter' )


# csv inputs
trueDE_dsn = data_path / "trueDE.csv"

counts_low_dsn = data_path / counts_low_csv
counts_high_dsn = data_path / counts_high_csv

edgeR_DESeq2_low_dsn = data_path / edgeR_DESeq2_low_csv 
edgeR_DESeq2_high_dsn = data_path / edgeR_DESeq2_high_csv 
 
######################################################################################################################################

method_stat_list =  ['edgeR_F', 'DESeq2_stat', 'bin_dev_F']
line_color_list = ['red','blue','black']


df_trueDE = pd.read_csv ( trueDE_dsn )[['x']].rename ( columns={'x':'index'} ).set_index( 'index' )
df_trueDE ['DE']=1
print (  '\n\n\n df_trueDE \n\n', df_trueDE )


df_counts_low = pd.read_csv ( counts_low_dsn, index_col=0 )
print (  '\n\n\n df_counts_low \n\n', df_counts_low )
 
df_counts_high = pd.read_csv ( counts_high_dsn, index_col=0 )
print (  '\n\n\n df_counts_high \n\n', df_counts_high )


df_DE_low = pd.read_csv ( edgeR_DESeq2_low_dsn, index_col=1 ).drop ( columns=['Unnamed: 0'] )
print (  '\n\n\n df_DE_low \n\n', df_DE_low )

df_DE_high = pd.read_csv ( edgeR_DESeq2_high_dsn, index_col=1 ).drop ( columns=['Unnamed: 0'] )
print (  '\n\n\n df_DE_high \n\n', df_DE_high )

######################################################################################################################################   

df_cell_data = pd.DataFrame ( index = df_counts_high.columns, data = 80*[0] + 80*[1], columns=['Group'] ) 
df_cell_data['Null_model'] = 0
n_cells = df_cell_data.shape[0]
df_cell_data['Saturated_model'] = range ( n_cells) 
print (  '\n\n\n df_cell_data: \n\n', df_cell_data )


print (  'calculate changes in binomial deviance for Lowly Expressed genes' ) 
df_changes_in_binomial_deviance_low,  df_NLL_low = calculate_changes_in_binomial_deviance ( df_counts_low )

print (  'calculate changes in binomial deviance for Highly Expressed genes' ) 
df_changes_in_binomial_deviance_high,  df_NLL_high = calculate_changes_in_binomial_deviance ( df_counts_high )
 
######################################################################################################################################   

df_compare = pd.concat ( [ df_changes_in_binomial_deviance_low, df_DE_low ], axis=1, sort=False )
df_AUC = plots ( 'Lowly'  )
df_AUC.to_excel(writer, sheet_name='Lowly expressed')

df_compare = pd.concat ( [ df_changes_in_binomial_deviance_high, df_DE_high ], axis=1, sort=False )
df_AUC = plots ( 'Highly'  )
df_AUC.to_excel(writer, sheet_name='Highly expressed')
  


pdf_pages.close()  
  
writer.save()






