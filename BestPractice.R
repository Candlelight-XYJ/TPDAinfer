#################################################################
## Using advanced TPDA algorithm to construct Bayesian Network ##
#################################################################
rm(list = ls())
setwd("E:/GitHub/TPDA/TPDAinfer/")

########################
## library R packages ##
########################
library(igraph)
library(graph)
library(bnlearn)
library(TPDAinfer)

##(1) data preprocessing
options(stringsAsFactors = F)
expMatrix <- read.csv("data/bayesinput.csv") # with header
expMatrix <- TPDAinfer::dfProcess(expMatrix)
##(2) get primary mutual information
gene2MI <- TPDAinfer::getMI(expMatrix, weight=0.01)
head(gene2MI)
##(3) start TPDA
weight1=0.01
weight2=0.01
weight3=0.01
filePath="data/bayesinput.csv"
TPDAresult <- TPDAinfer::TPDA_algorithm(filePath,
                         weight1,
                         weight2,
                         weight3,
                         gene2MI)
##(4) construct Bayesian network in R
BN <- TPDAinfer::convertBN(TPDAresult)
##(5) Parameter Learning
set.seed(1000) # setting seeds
## data discretize
expData <- bnlearn::discretize(data.frame(t(expMatrix)), method ='quantile', breaks=2 )
rownames(expData) <- colnames(expMatrix)
fitted_time<-system.time(fitted<-bnlearn::bn.fit(BN,data = expData,method='bayes')) # bayes only for discrete data
fitted_time # show learning times

#############################
##    Bayes Inferrence     ##
#############################
## Every gene status
nodeStatus <- data.frame()
for(i in 1:ncol(expData)){
  print(i)
  #i=1
  status <- as.character(unique(expData[,i]))
  nodeStatus <- rbind(nodeStatus, c(colnames(expData)[i],status))
}
head(nodeStatus)

## (2)bayes inferrence
## Calculating the conditional probability of Gene 1 and Gene 2 for phenotype 1
cpquery(fitted,
        event = (pheno1=='(0.015,2.97]'),
        evidence =((gene1=='[-3,-0.00341]')&(gene2=='(-0.0443,3]')))

## [1] 0.4894299
## [1] 0.4939955

