
import pandas as pd
import numpy  as np


import sys
sys.path.append("D:/scRNA-seq_DE/DE_scRNAseq")
 
from  utilities import *



from random import shuffle




#########################################################################################
#########################################################################################


def  compute_binomial_deviance__genuine_and_randomized_counts ( df_UMI_counts, n_randomizations=11, pctl_list =  [.25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ] ):
  
  df_genuine = compute_binomial_deviance ( df_UMI_counts, 'genuine' )

  UMI_counts_index = df_UMI_counts.index
  df_deviance_randomized_list = []

  for rep in range ( n_randomizations ):
    print  ( 'randomization: ', rep )
  
    df_UMI_counts_randomized_list= []

    for col in df_UMI_counts.columns.values.tolist() :
      df_col_list =  df_UMI_counts [[ col ]].values.tolist()
      np.random.shuffle(df_col_list)
      df_col_randomized = pd.DataFrame ( index=UMI_counts_index, data = df_col_list, columns = [ col ] )
  
      df_UMI_counts_randomized_list.append ( df_col_randomized )

    df_UMI_counts_randomized = pd.concat ( df_UMI_counts_randomized_list, axis=1, sort=True )

    df_deviance_randomized = compute_binomial_deviance ( df_UMI_counts_randomized, rep )  
  
    df_deviance_randomized_list.append ( df_deviance_randomized )
    df_deviance_randomized_multis = pd.concat ( df_deviance_randomized_list,  axis=1, sort=True )


  
  df_percentiles_randomized = df_deviance_randomized_multis.describe( percentiles=pctl_list )   
  sel_index_list = df_percentiles_randomized.transpose().drop ( columns = [ 'count', 'mean', 'std', 'min' ] ).columns.values.tolist()  
  
  ser_median_percentiles = df_percentiles_randomized.median( axis = 1 )  
  
  
  df_tuple_list = []  
  for percentile in sel_index_list :
    threshold = ser_median_percentiles.loc [ percentile ]
    sel_gene_count = np.sum ( ( df_genuine ['genuine' ] >=  threshold )  )
    df_tuple_list.append ( ( percentile, threshold, sel_gene_count ) )
  
  df_summary_table = pd.DataFrame ( data = df_tuple_list, columns= ['percentile rank', 'median percentile - randomized data', '# genes selected' ] ).sort_index ( ascending=False ).set_index ( ['percentile rank'] )
  
  return  [ df_genuine, df_deviance_randomized_multis, df_summary_table ]
    
  
  
  
def compute_binomial_deviance ( df_UMI_counts, column_name ):  

  arr_gen_counts = df_UMI_counts.values 
  mat_gen_counts = np.matrix ( arr_gen_counts)

  n_vector = mat_gen_counts.sum(axis=0)
  p_vector = mat_gen_counts.sum(axis=1)

  pi_hat_vector = p_vector/ p_vector.sum()
  n_pi_hat = np.array ( pi_hat_vector * n_vector )
  term1 = arr_gen_counts * np.log ( 1e-20 + np.nan_to_num( arr_gen_counts / n_pi_hat ) )

  arr_ni_m_gen_counts = np.array( n_vector ) - arr_gen_counts
  n_1_m_pi_hat =  np.array( n_vector ) - n_pi_hat
  term2 = arr_ni_m_gen_counts * np.log ( 1e-20 +  arr_ni_m_gen_counts / n_1_m_pi_hat )

  rhs = 2 * ( term1 + term2 )


  df_deviance = pd.DataFrame( index = df_UMI_counts.index, data = rhs.sum( axis=1 ), columns=[column_name] ) 

  return df_deviance
 
 
################################################################################################################################
  
def filter_UMI_counts ( df_UMI_counts, df_binomial_deviance, number_selected_genes ):


  df_binomial_deviance['rank'] = df_binomial_deviance['genuine'].rank ( ascending=False )
  df_bd_select = df_binomial_deviance.loc [ df_binomial_deviance['rank'] <= number_selected_genes ]
  df_filtered_UMI_counts = df_UMI_counts.loc [ df_bd_select.index ] 
  
  return  [ df_filtered_UMI_counts ]
  
#########################################################################################  

def  MLE_cell_cluster_models  ( df_UMI_counts, df_cell_clusters ):

  n_cells = df_cell_clusters.shape[0] 
  n_genes = df_UMI_counts.shape[0]

  df_cell_clusters[1] = 0
  df_cell_clusters[ n_cells ] = range ( n_cells )


  clusterings_list = df_cell_clusters.columns.values.tolist()
  clusterings_list.sort(    )
  dict_MLE_results = {}
  NLL_list = []
  

  for clustering in  clusterings_list:
   
    df_cell_cluster_gene_totals_list = []
  
    for cluster in range ( clustering ):
      cluster_cell_list = df_cell_clusters[[clustering]].loc [ df_cell_clusters[clustering] ==cluster ].index.values.tolist()
      df_UMI_counts_cluster = df_UMI_counts [ cluster_cell_list ]
      df_cell_cluster_gene_totals_list.append ( df_UMI_counts_cluster.sum ( axis=1 ).to_frame ( name = cluster ) )  
    
    df_cell_cluster_gene_total_counts = pd.concat ( df_cell_cluster_gene_totals_list, axis=1, sort=False ) 

    ser_cluster_total_counts = df_cell_cluster_gene_total_counts.sum() 
    df_pi_hat = df_cell_cluster_gene_total_counts. divide ( ser_cluster_total_counts )
  
    df_product_count_log_pi_hat = df_cell_cluster_gene_total_counts.multiply ( np.log ( 1.e-20 + df_pi_hat ) )
    NLL = - df_product_count_log_pi_hat.sum().sum()  
    NLL_list.append ( NLL ) 
 
    dict_cluster_results =  { 'cell_cluster_gene_totals': df_cell_cluster_gene_total_counts, 'pi_hat':df_pi_hat, 'NLL': NLL  }
    dict_MLE_results [ clustering ] =  dict_cluster_results   

	

  df_NLL = pd.DataFrame ( index = clusterings_list )
  df_NLL ['NLL'] = NLL_list
  df_NLL['LR statistic'] = np.nan

  NLL_null = df_NLL.at [ 1, 'NLL']
  for cluster in  clusterings_list [1:]:
    NLL = df_NLL.at [ cluster, 'NLL' ] 
    df_NLL.at [ cluster, 'LR statistic'] = 2 * ( NLL_null - NLL )
 
  return [ dict_MLE_results, df_NLL ] 

#########################################################################################  
 
  