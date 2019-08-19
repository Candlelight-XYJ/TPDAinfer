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
expMatrix <- read.csv("../testData/bayesinput.csv") # with header
expMatrix <- dfProcess(expMatrix)
##(2) get primary mutual information
gene2MI <- getMI(expMatrix, weight=0.01)
head(gene2MI)
##(3) start TPDA
weight1=0.01
weight2=0.01
weight3=0.01
filePath="../testData/bayesinput.csv"
TPDAresult <- TPDA_algorithm(filePath,
                         weight1,
                         weight2,
                         weight3,
                         gene2MI)
##(4) construct Bayesian network in R
BN <- convertBN(TPDAresult)
##(5) Parameter Learning
set.seed(1000)## setting seeds
expData <- data.frame(t(expMatrix))
fitted_time<-system.time(fitted<-bn.fit(BN,data = expData,method='mle')) #bayes only for discrete data
fitted_time ## show learning times

#############################
##    Bayes Inferrence     ##
#############################
## For Inferrence ,  data must be discrete
expData <-

## (1)preprocess quantile data
lisan<-data.frame()
for(i in 1:ncol(expData))
{
  xx<-unique(mydata[,i])
  lisan<-rbind(lisan,t(as.data.frame(xx)))
}
add_rowname<-read.csv(input_path)
all_rowname<-add_rowname[,1]                    ##add gene names
##no_repeat<-read.csv(output_lisan_path)
no_repeat<-lisan
no_repeat$X<-add_rowname[,1]                    ##combine gene name and quantile intervals
head(no_repeat)

## (2)bayes infer
## Calculating the reasoning probability of Gene 1 and Gene 2 for phenotype 1
cpquery(fitted,
        event = (pheno1=='(0.015,2.97]'),
        evidence =((gene1=='[-3,-0.00341]')&(gene2=='(-0.0443,3]')))

## [1] 0.5070332


