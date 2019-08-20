## A bayesian network workflow to mine gene-gene and gene-pheno relationships


## Publication & How to Cite
+ The TPDApipe publication can be found in: http://www.ncbi.nlm.nih.gov/pubmed/xxxxxx

+ If you used this pipeline for your analysis, please cite: XXX, DOI:XXXXXXX

+ Thanks in advance!

## Files including
+ An R package `TPDAinfer`
+ Bestpractice.R

## Best Practice

#### 1) Installation
+ Install the latest development version from GitHub with
```r
devtools::install_github("Candlelight-XYJ/TPDA")
```

#### 2) Data preprocessing
+ Using `dfProcess` function to preprocess expression data
```r
rm(list = ls())
library(igraph)
library(graph)
library(bnlearn)
library(TPDAinfer)
options(stringsAsFactors = F)
expMatrix <- read.csv("data/bayesinput.csv") # with header
expMatrix <- TPDAinfer::dfProcess(expMatrix,discrete=F)
```
#### 3) Computating mutual information
+ Preparing primary mutual information for next structure learning
```r
gene2MI <- TPDAinfer::getMI(expMatrix, weight=0.01)
head(gene2MI)
```

#### 4) Using advanced TPDA algorithms to construct Bayesian network
```r
weight1=0.01
weight2=0.01
weight3=0.01
filePath="data/bayesinput.csv"
TPDAresult <- TPDAinfer::TPDA_algorithm(filePath,
                                         weight1,
                                         weight2,
                                         weight3,
                                         gene2MI)
```

#### 5) Convert data.frame to graph format
+ Converting TPDA network result to bnlearn graph format in R 
```r
BN <- TPDAinfer::convertBN(TPDAresult)
```

#### 6) Parameter learning
```r
set.seed(1000) # setting seeds
## data discretize
expData <- bnlearn::discretize(data.frame(t(expMatrix)), method ='quantile', breaks=2 )
rownames(expData) <- colnames(expMatrix)
fitted_time<-system.time(fitted<-bnlearn::bn.fit(BN,data = expData,method='bayes')) # bayes only for discrete data
fitted_time # show learning times
```

#### 7) Bayes inferrence
+ Every gene status
```r
nodeStatus <- data.frame()
for(i in 1:ncol(expData)){
  print(i)
  #i=1
  status <- as.character(unique(expData[,i]))
  nodeStatus <- rbind(nodeStatus, c(colnames(expData)[i],status))
}
head(nodeStatus)
```

+ Bayes inferrence
+ Calculating the conditional probability of Gene 1 and Gene 2 for phenotype 1
```r
cpquery(fitted,
        event = (pheno1=='(0.015,2.97]'),
        evidence =((gene1=='[-3,-0.00341]')&(gene2=='(-0.0443,3]')))

## [1] 0.4894299
## [1] 0.4939955
```

