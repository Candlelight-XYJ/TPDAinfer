#########################################
## Bayesian Network Parameter Learning ##
#########################################

##########################
##     files‘ paths     ##
##########################
input_path <- "E:\\GitHub\\TPDA\\testData\\bayesinput.csv"
TPDA_structure <- "E:\\GitHub\\TPDA\\outputData\\2_TDPA_structure.csv"
############################
##    Data Preprocess     ##
############################
dfProcess <- function(input_path){
  dataframe<-read.csv(input_path,header = F,stringsAsFactors = FALSE)[-1,]
  row.names(dataframe)<-dataframe[,1]
  dataframe<-dataframe[-1]
  mydata<-data.frame(dataframe)  
  tq <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
  rownames(tq) <- rownames(mydata)
  colnames(tq) <- colnames(mydata)
  mydata <- tq
  mydata<-data.frame(t(mydata))
  tpN <- discretize(mydata,method ='quantile',breaks=2 ) #离散化
  rownames(tpN) <- rownames(mydata)
  colnames(tpN) <- colnames(mydata)
  mydata <- data.frame((tpN))  
  return(mydata)
}
## 
zp <- dfProcess(input_path)

#######################################################################
##  Using TPDA structure and graphNEL() to construct Network in R    ##
#######################################################################

library(graph)
result<-read.csv(TPDA_structure,stringsAsFactors = FALSE) ## input TPDA algorithm result
v<-union(result[,1],result[,2])
edL=vector("list",length=length(v)) ## length(v) means nodes num
names(edL)<-v
for(i in 1:length(v))
{   edge_num<-0
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

gRbn<-as.bn(gR)  ##convert gR to bnlearn object


###############################
##    Parameter Learning     ##
###############################
set.seed(1000)##设置种子点
fitted_time<-system.time(fitted<-bn.fit(gRbn,data = mydata,method='bayes')) 理
fitted_time ## show learning times


