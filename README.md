This set of programs analyzes differential expression for single cell RNA-Seq data, represented by UMI counts.
<br />
<br />
The approach closely follows the work of Townes et al. [1] by using
- multinomial models
- binomial deviance to filter genes for clustering
<br />



Consider a cell clustering with K clusters
- For each cluster k, the multinomial model gives maximum likelihood estimates of the relative abundance of each gene.  
- Denote the estimate for gene j as  pi_hat(k, j)
<br />
For each gene j, the range of values
<br />
<br />
r(j) =  max ( pi_hat(k,j)  ) -  min ( pi_hat(k,j)  )  -- max and min taken over k
<br />
<br />
measures differential expression and is correlated with binomial deviance.
<br />
<br />
For details, please refer to 
- example_programs.docx
- the folder example_programs_Zhengmix8eq: a 6-program stream to analyze the Zhengmix8eq data set
- the folder example_programs_simulated: a 4-program stream to create and analyze a simulated data set with 2 clusters
- functions that perform the analysis, and documentation:
  - FUNCTIONS_DE_scRNAseq.py
  - utilities.py
  - function documentation.docx
<br />
<br />
Reference
<br />
1. Townes F W, Hicks S C, Aryee M J  et al.: Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model. Genome Biol 20, 295 (2019). https://doi.org/10.118s13059-019-1861-6 
