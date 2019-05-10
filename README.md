## What is TPDApipe


## Publication & How to Cite
+ The TPDApipe publication can be found in: http://www.ncbi.nlm.nih.gov/pubmed/xxxxxx

+ If you used this pipeline for your analysis, please cite: YuJia Xiang,JianXiao Liu (2019). TPDAPipe: a bayesian network workflow to mining gene-gene relationships. Bioinformatics, DOI:XXXXXXX

+ Thanks in advance!

## Files including
+ 0-compute_mutual_info.py
+ 1-TPDA_construct_BN.R
+ 2-parameter_learning_and_BN_infer.R


## Best Practice

#### 1) Computating mutual information
By using our workflow,the first thing to do is to computate mutual information of gene-gene or gene-pheno .

#### 2) Using advanced TPDA algorithms to construct Bayesian network
We provide the source code of TPDA algorithm (Script `1-TPDA_construct_BN.R`),to help user construct Bayesian network .
Open the script ,you need provide 6 parameters
+ weight_1, weight_2, weight_3

This is the threshold for user to set,different weight will construct different BN structure. The default weight is 0.01, We suggest user use proper weight for your special data.

+ input_path
This is the file path, you input your gene expression data file

+ gene2MI_path
This is a Output path, Outputting primary MI values

+ network_structure_path
This is a Output path, Outputting the BN result of TPDA

#### 3) Parameter learning



#### 4) Bayesian inferrence
