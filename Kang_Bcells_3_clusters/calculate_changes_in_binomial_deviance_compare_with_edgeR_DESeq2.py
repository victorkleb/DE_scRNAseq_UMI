

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





pd.options.display.width = 180
pd.set_option('display.max_columns', 30)

#########################################################################################

def abs_true_DE ( df_row_data, n_groups ):
  print ( 'df_row_data: \n\n', df_row_data )
  group_list = list ( range ( n_groups ) )   
  DEFacGroup_col_list = [ 'DEFacGroup' + str ( i+1) for i in group_list ]
  rn_dict = dict ( zip ( DEFacGroup_col_list, group_list ) )    
  df_DEFacGroup_cols = df_row_data [ DEFacGroup_col_list ].rename ( columns= rn_dict )

  df_true_DE_list = []  
  col_count = 0  
  for i  in range ( 1, n_groups ):
    for j in range ( i ):
      df_temp = df_DEFacGroup_cols [[ i, j ]]
      df_temp [ col_count ] = np.abs ( np.log2 ( df_temp[i] ) - np.log2 ( df_temp[j] ) )
      df_true_DE_list.append ( df_temp[ col_count ] )
      col_count += 1	  
   		
  df_true_DE_all = pd.concat ( df_true_DE_list, axis=1, sort=False )
  ser_abs_true_DE = df_true_DE_all.max ( axis=1 ) 
  
  return  ser_abs_true_DE



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
  
  deg_freedeom_num = n_groups - 1
  deg_freedeom_denom = n_cells - n_groups  
  df_changes['bin_dev_F'] = ( df_changes['cluster_to_null'] / deg_freedeom_num ) / ( df_changes['saturated_to_cluster'] / deg_freedeom_denom )

  df_changes['log10_bin_dev_F'] = np.log10 ( df_changes['bin_dev_F'] )
  
  df_changes['bin_dev_p_value'] =  df_changes['bin_dev_F'].apply ( lambda F: 1 - f.cdf( F, deg_freedeom_num,  deg_freedeom_denom ) )

  reject, padj_array = fdrcorrection( df_changes['bin_dev_p_value'].values, alpha=0.05  )
  df_changes['bin_dev_p_value_adj'] = padj_array  
  print (  '\n\n\n df_changes: \n\n', df_changes  ) 

  return ( df_changes, df_NLL )


  

def scatter_plot ( df_plot, x, y, x_label, y_label ):

  fig, ax1 = plt.subplots()  
  titl = title1 +  title2 + title3
  plt.title( titl, fontsize=5.5 )

  ax1.scatter ( df_plot[x], df_plot[y],  c='black', s=2 ) 
  
  ax1.set_xlabel ( x_label, fontsize=6.1 )	
  ax1.set_ylabel ( y_label, fontsize=6.1 )	
  ax1.tick_params(labelsize=5.2) 
 
  pdf_pages.savefig( fig, transparent=True )
  

  

def ROC_plots ( DE_threshold ):
  global title1, title2
  
  df_compare_methods_notNA['DE'] = ( df_compare_methods_notNA['abs_true_DE'] > DE_threshold ).astype(int)
  number_DE = df_compare_methods_notNA['DE'].sum()  
  
  dict_method_ROCs = {}
  for stat,color in zip ( method_stat_list, line_color_list ):
    fpr, tpr, thresholds = metrics.roc_curve(  df_compare_methods_notNA['DE'].values,  np.abs ( df_compare_methods_notNA[stat].values ) ) 
    df_ROC = pd.DataFrame ( data = zip ( fpr,tpr ), columns = ['fpr','tpr'] )    
    AUC = metrics.roc_auc_score(  df_compare_methods_notNA['DE'].values,  np.abs ( df_compare_methods_notNA[stat].values ) ) 
    dict_method_ROCs[ stat ] = { 'df_ROC':df_ROC, 'AUC': AUC, 'color':color }       
 
    
  title3 = '\n ROC curves for edgeR, DESeq2, and bin_dev_F'
  title4 = '\n define genes as Differentially Expressed if  max ( abs ( log2 ( DEFacGroupI / DEFacGroupJ ) ) ) > '  + '{0:.2f}'.format( DE_threshold )  
  title5 = '\n  number Differentially Expressed: ' +    '{0:.0f}'.format( number_DE )
 
  
  fig, ax1 = plt.subplots()
  titl = title1 +  title2 + title3 + title4 + title5
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
  df_AUC['DE_threshold'] = DE_threshold
  
  return  df_AUC   
 
#########################################################################################
 
work_folder =  "D:/scRNA-seq_DE/simulation_Kang_Bcells_3_clusters"
data_path = Path ( work_folder )

AUC_DE_simulations_dsn = data_path / "AUC_DE_simulations.pkl"


##  PDF  base dsn
plot_base = "compare_methods_"

# Excel output base
Excel_base = "compare_methods_"
 
 
# csv output base
compare_base = "compare_bin_dev_F_with_edgeR_DESeq2_" 


# csv input base dsns 
splatter_counts_base  = "splatter_counts_" 
splatter_rows_base = "splatter_rows_"
splatter_columns_base = "splatter_columns_"
edgeR_DESeq2_base     = "splatter_DE_"   

######################################################################################################################################

number_of_simulations = 10

 

method_stat_list =  ['edgeR_F', 'DESeq2_stat', 'bin_dev_F']
line_color_list = ['red','blue','black']

df_AUC_list = []

for simulation in  range(1, 1 + number_of_simulations ):
  print (  '\n\n simulation: ', simulation )
  
  str_sim = str( simulation )

  row_data_csv = splatter_rows_base + str_sim + ".csv" 
  column_data_csv = splatter_columns_base + str_sim + ".csv"
  counts_csv = splatter_counts_base + str_sim + ".csv"
  edgeR_DESeq2_csv = edgeR_DESeq2_base + str_sim + ".csv"
  compare_csv = compare_base + str_sim + ".csv"  
  plot_pdf = plot_base + str_sim + ".pdf"

  row_data_dsn = data_path / row_data_csv  
  column_data_dsn = data_path / column_data_csv
  counts_dsn = data_path /  counts_csv
  edgeR_DESeq2_dsn = data_path / edgeR_DESeq2_csv
  compare_dsn = data_path / compare_csv  
  plot_dsn = data_path / plot_pdf
  pdf_pages = PdfPages( plot_dsn )  

  #####  groups are numbered 1,2 ;  must change to 0,1  for python and MLE computations
  df_cell_data = pd.read_csv ( column_data_dsn, index_col=0 ) [['Group']].rename ( columns={'Group':'splatter_group'} )
  df_cell_data['Group'] = df_cell_data['splatter_group'].str.slice(5, 6) .astype(int) - 1
  df_cell_data['Null_model'] = 0
  n_cells = df_cell_data.shape[0]
  df_cell_data['Saturated_model'] = range ( n_cells) 
  print (  '\n\n\n df_cell_data: \n\n', df_cell_data )

  n_groups = 1 +  df_cell_data['Group'].max()
  print (  '\n\n\n n_groups: \n\n', n_groups )  
  

  df_row_data = pd.read_csv ( row_data_dsn, index_col=0 ).drop ( columns=['Gene'] )
  df_row_data['abs_true_DE'] = abs_true_DE ( df_row_data, n_groups )       
  print ( '\n\n\n df_row_data: \n\n', df_row_data )  



  
  df_counts = pd.read_csv ( counts_dsn, index_col=0 ) 
  print (  '\n\n\n df_counts: \n\n', df_counts )

  df_edgeR_DESeq2 = pd.read_csv ( edgeR_DESeq2_dsn, index_col=1 ).drop ( columns=['Unnamed: 0'] )
  df_edgeR_DESeq2['log10_edgeR_F'] = np.log10 ( df_edgeR_DESeq2['edgeR_F'] )  
  df_edgeR_DESeq2['log10_DESeq2_stat'] = np.log10 ( df_edgeR_DESeq2['DESeq2_stat'] )  
  print (  '\n\n\n df_edgeR_DESeq2 \n\n', df_edgeR_DESeq2 )


  print (  '\n\n\n\n calculate changes in binomial deviance' ) 
  df_changes_in_binomial_deviance,  df_NLL = calculate_changes_in_binomial_deviance ( df_counts )

  
  df_compare = pd.concat ( [ df_row_data , df_edgeR_DESeq2, df_changes_in_binomial_deviance ], axis=1, sort=False ) 
  print (  '\n\n\n df_compare: \n\n', df_compare )   
   
 
  df_compare_methods_notNA = df_compare.loc [ df_compare['bin_dev_F'].notna() ]

  df_compare.to_csv ( compare_dsn )    
  
  
  n_genes = df_compare.shape[0]  
       
  
  title1 = "Splatter simulated data with " +  str( n_groups )  + " clusters - data set " + str_sim 
  title2 = "\n parameters derived from " + str ( n_genes ) + " most highly expressed genes - Kang Lupus B cell data"
  
  title3 = '\n compare DESeq2 statistic to true DE parameter'
  scatter_plot ( df_compare_methods_notNA,  'abs_true_DE', 'log10_DESeq2_stat', 'true Differential Expression parameter', 'Log10 ( DESeq2 statistic )' )

  title3 = '\n compare edgeR F-statistic to true DE parameter'   
  scatter_plot ( df_compare_methods_notNA ,  'abs_true_DE', 'log10_edgeR_F', 'true Differential Expression parameter', 'Log10 ( edgeR F-statistic )' )

  title3 = '\n compare binomial deviance F-statistic to true DE parameter'   
  scatter_plot ( df_compare_methods_notNA ,  'abs_true_DE', 'log10_bin_dev_F', 'true Differential Expression parameter', 'Log10 ( binomial deviance F statistic )' )
    
  
  
  df_AUC_DE_list = []  
  for  DE_threshold in [ 1, .8, .6, .4, .2, .1 ]:
    df_AUC_DE_list. append ( ROC_plots ( DE_threshold ) )
  
  df_AUC_DE = pd.concat ( df_AUC_DE_list, sort=False )
  df_AUC_DE['simulation'] = simulation
  print ( '\n\n\n df_AUC_DE: \n\n', df_AUC_DE )

  df_AUC_list.append ( df_AUC_DE )  

  pdf_pages.close()  
  
  
       
  
df_AUC = pd.concat ( df_AUC_list, sort=False )
print ( '\n\n\n df_AUC_: \n\n', df_AUC )



df_AUC.to_pickle ( AUC_DE_simulations_dsn )

 
  
 





