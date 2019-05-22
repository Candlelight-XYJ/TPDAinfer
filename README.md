## A bayesian network workflow to mine gene-gene and gene-pheno relationships


## Publication & How to Cite
+ The TPDApipe publication can be found in: http://www.ncbi.nlm.nih.gov/pubmed/xxxxxx

+ If you used this pipeline for your analysis, please cite: XXX, DOI:XXXXXXX

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
+ **weight_1, weight_2, weight_3**

This is the threshold for user to set,different weight will construct different BN structure. The default weight is 0.01, We suggest user use proper weight for your special data.

+ **input_path**

This is the file path, you input your gene expression data file

+ **gene2MI_path**

This is a Output path, Outputting primary MI values

+ **network_structure_path**

This is a Output path, Outputting the BN result of TPDA

After providing file paths,you can do as follows:

+ **`First`** load Preprocess Functions
  + dfProcess <- function(input_path)
  + getDoubleHang_value<-function(mylist1,mylist2)
  + getMI <- function(df, weight)
  + getCMI <- function(gxi, gyi,cutset,gen)

+ **`Second`** load Three Phase Development Algorithm Function

TPDA_algorithm <- function(weight_1,weight_2,weight_3,
                           input_path,
                           gene2MI_path)



#### 3) Parameter learning
We provide Script `2-parameter_learning_and_BN_infer.R` for user to learn parameters for their Bayesian network. 
+ **`First`** input your files path
```r
input_path <- "E:\\GitHub\\TPDA\\testData\\bayesinput.csv"
TPDA_structure <- "E:\\GitHub\\TPDA\\outputData\\2_TDPA_structure.csv"
```
+ **`Second`**  using  Data Preprocess function to format input data
```r
mydata <- dfProcess(input_path)
```
+ **`Third`**  Using TPDA structure result and graphNEL() to construct Network in R
```r
library(graph)
result<-read.csv(TPDA_structure,stringsAsFactors = FALSE) ## input TPDA algorithm result
v<-union(result[,1],result[,2])
edL=vector("list",length=length(v)) ## length(v) is the num of gene nodes
names(edL)<-v
for(i in 1:length(v))
{ edge_num<-0
nodename<-names(edL[i])
node_edge<-c()
for(j in 1:nrow(result))
{ 
  if(as.character(result[j,1])==nodename) 
  {       
    son_node<-as.character(result[j,2])
    node_edge<-c(node_edge,son_node)
  } 
}
edL[[i]]<-list(edges=node_edge)
}
gR <- graphNEL(nodes=v, edgeL=edL,edgemode = "directed")
gRbn<-as.bn(gR) ##convert gR to bnlearn object

```

+ **`Fourth`** Parameter Learning
Using function `bn.fit` in R package bnlearn, to learn you BN parameters
```r
set.seed(1000)## setting seeds
fitted_time<-system.time(fitted<-bn.fit(gRbn,data = mydata,method='bayes')) 
fitted_time ## show learning times
```

#### 4) Bayesian inferrence
In script 
`2-parameter_learning_and_BN_infer.R`
