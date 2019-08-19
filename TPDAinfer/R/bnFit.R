########################################
##                                    ##
## code last update: August 19, 2019  ##
##                                    ##
########################################

##' Converting TPDA structure result to bnlearn network format
##' @param TPDAstruct
##' @return bayesian network structure
##' @import graph, bnlearn
##' @description Using TPDA structure result and graphNEL() to construct Network in R
convertBN <- function(TPDAstruct){
  v<-union(TPDAstruct[,1],TPDAstruct[,2])
  edL=vector("list",length=length(v)) ## length(v) is the num of gene nodes
  names(edL)<-v
  for(i in 1:length(v))
  {   edge_num<-0
  nodename<-names(edL[i])
  node_edge<-c()
  for(j in 1:nrow(TPDAstruct))
  {
    if(as.character(TPDAstruct[j,1])==nodename)
    {
      son_node<-as.character(TPDAstruct[j,2])
      node_edge<-c(node_edge,son_node)
    }
  }
  edL[[i]]<-list(edges=node_edge)
  }
  gR <- graphNEL(nodes=v, edgeL=edL,edgemode = "directed")
  gRbn<-as.bn(gR)  ##convert gR to bnlearn object
  return(gRbn)
}



