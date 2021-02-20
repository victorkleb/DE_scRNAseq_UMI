## DE_scRNAseq_UMI

The ratio of **changes in binomial deviance** between apprpriate multinomial models identifies Differentially Expressed genes.  This approach draws on the work of Townes et al. [1]  in which binomial deviance is proposed to identify informative genes for clustering. 
<br /><br />
Given 
-	**A**  a g by c array of counts, with g genes and c cells
-	in which cells are grouped in K clusters   C_1 – C_K
<br />

Define three multinomial models for the relative abundance of each gene
<br />
-	Null model
-	Saturated model
-	Cluster model
<br />

The maximum likelihood estimates of the model parameters are

-	Null model: for each gene <br />
pi_hat_Null ( g ) = sum ( **A**( g , . ) ) / sum ( **A**( . , . ) )

-	Saturated model:  for each gene and cell  <br />
pi_hat_Sat ( g, c )  = **A**( g, c ) / sum ( **A** ( . , c ) ) 

-	Cluster model:  for each gene and cluster  <br />
pi_hat_Cluster ( g, k ) = sum( ( **A**( g, c in C_k  ) ) / sum( ( **A**( . , c in C_k  ) )
<br />

Similar to the calculation of binomial deviance of the Null model in the section “Feature selection using deviance” in [1], the change in binomial deviance is calculated 
-	between the Cluster and Null models
-	between the Saturated and Cluster models (that is, the binomial deviance of the Cluster model)
<br />

The distributions of these are approximately chi-square.  Hence, their ratio is approximately an F statistic, which is proposed as a measure of Differential Expression.  For convenience, this method is referred to as bin_dev_F.
<br />
<br />
This posting includes examples with simulated data.  For all examples that have been studied, bin_dev_F is comparable to **edgeR** and **DESeq2**.
<br />
<br />

**Example 1:  Mou et al.**

The study by Mou et al. [2] includes simulated data [3] that may be downloaded from GitHub.  

The two simulated data sets differ in the fold change parameter for the degree of DE.  Here, results are given for the data set with the smaller fold change parameter.     

The authors identify highly and lowly expressed genes.  The two sets were analyzed separately.  Evaluating performance with AUC – area under the receiver operating characteristic (ROC) curve – shows that bin_dev_F is comparable to **edgeR** and **DESeq2**.
<br />
<br />	
The posted folder Tianmou_1 contains
- an R program to 
  - read the Rdata workspace downloaded from GitHub
  - output counts and metadata to CSV files
  - perform DE with **edgeR** and **DESeq2**, with CSV output
- a python program to
  - perform DE with bin_dev_F
  - compare the three methods’ results
- a PDF file containing plots for the sets of highly and lowly expressed genes
  - scatter plots comparing the 3 DE methods’ statistics
  - ROC curves
<br />
<br />

**Example 2:  Simulated data with parameters derived from Kang Lupus data**

The Bioconductor web page for the **muscat** package [4] includes a vignette [5].   

Following this vignette, highly expressed genes for B cells were selected to provide parameters for the **splatter** simulation package [6,7,8].  Ten sets of simulated data, each with two clusters, were generated with **splatter**.  
<br />
Each simulated data set was analyzed with **edgeR**, **DESeq2**, and bin_dev_F.
<br />
<br />
**Splatter** does not produce a dichotomous DE indicator.  Instead, for each gene it outputs values DEFac[Group]  described in the documentation as “the differential expression factor for each gene in a particular group.”  In the two-cluster setting, true DE parameters can be obtained from these as 
<br />
<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; DE_parameter = log2 (  DEFacGroup2 /  DEFacGroup1 )

<br />
The results of each DE detection method were compared to these true parameters in scatter plots. 
<br />
<br />
As with the Mou et al. data, DE methods were compared by calculating the area under the ROC curve.  Since this requires a dichotomous classification for each gene (that is: is the gene DE or not) the true DE parameters were dichotomized by thresholding.  Given a threshold, define
<br />
<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; DE_classification = (  abs ( DE_parameter ) > threshold ) 
<br />
<br />
AUC was calculated for 50 conditions
<br />

- 10 simulated data sets
- 5 thresholds
<br />

For each condition, bin_dev_F is compared with **edgeR** and **DESeq2**  by calculating the differences
<br />

- AUC ( bin_dev_F ) – AUC ( **edgeR** )
- AUC ( bin_dev_F ) – AUC ( **DESeq2** )
<br />
For the 50 conditions, bin_dev_F performance is
<br />

- almost as good as **edgeR**
- better than **DESeq2**
<br />
<br />	
  
The posted folder Kang_Bcells  contains
- an R program that
  - follows the muscat Bioconductor vignette [5]
  - prepares simulated data with **splatter**
  - outputs counts and metadata to CSV files
  - performs DE with **edgeR** and **DESeq2**, with CSV output
- a python program to
  - perform DE with bin_dev_F
  - compare the three methods’ results
- two PDF files
  - scatter plots and ROC curves for a single simulated data set
  - violin plots of the difference of AUC between bin_dev_F and both **edgeR** and **DESeq2** for the 50 data conditions analyzed
<br />
<br /> 	 

**References**
<br />
 
1. Townes F, Hicks S., Aryee M, et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model. Genome Biol 20, 295 (2019). https://doi.org/10.1186/s13059-019-1861-6

2. Mou T, Deng W, Gu F, Pawitan Y, et al. Reproducibility of Methods to Detect Differentially Expressed Genes from Single-Cell RNA Sequencing, Frontiers in Genetics, 10, 1331 (2020).
https://www.frontiersin.org/article/10.3389/fgene.2019.01331     

3. https://github.com/Tianmou/scRNAseq-DE-comparison

4. Crowell H, Germain P, Soneson C, Sonrel A, Robinson M.  muscat: Multi-sample multi-group scRNA-seq data analysis tools. R package version 1.4.0, https://github.com/HelenaLC/muscat.

5. https://bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html

6. Zappia L, Phipson B, Oshlack A.  Splatter: simulation of single-cell RNA sequencing data. Genome Biol 18, 174 (2017). https://doi.org/10.1186/s13059-017-1305-0

7. http://www.bioconductor.org/packages/release/bioc/html/splatter.html

8. http://www.bioconductor.org/packages/release/bioc/vignettes/splatter/inst/doc/splatter.html
 
