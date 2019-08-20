########################################
##                                    ##
## code last update: August 19, 2019  ##
##                                    ##
########################################

##' Three Phase Development Algorithm
##'
##' @param filePath a file path of gene expression data
##' @param weight1,weight2,weight3 weight values, mutual info threshold for  phase1, phase2, phase3
##' @param gene2MI a data frame containing gene-gene mutual information
##' @return bayesian network structure
##' @export
TPDA_algorithm <- function(filePath, weight1, weight2, weight3, gene2MI){
  expMatrix<-read.csv(filePath,header = FALSE,stringsAsFactors = FALSE)
  ck <- union(gene2MI$from,gene2MI$to)
  lname<-expMatrix[,1]
  expMatrix<-expMatrix[-1,]
  labe<-intersect(lname,ck)

  #covert n.x to gene.node
  Translat <- function(data_nx){
    retdata <- data_nx
    for(i in 1:length(data_nx)){
      retdata[i] <- labe[which( trans_labe == data_nx[i])]
    }
    return(retdata)
  }
  #covert gene object to n.x object
  trans_labe<-c()
  for(i in 1:length(labe)){ trans_labe<-c(trans_labe,paste("n.",i,sep = "")) }
  ct<-c()
  snum<-length(gene2MI$from)
  for(i in 1:snum){ ct<-c(ct,gene2MI[i,1],gene2MI[i,2]) }
  for(i in 1:length(ct)){
    #print(paste(i,ct[i]))
    ct[i]<-paste("n.",which(labe == ct[i]),sep = "")
  }

  ## Start TPDA
  ##=====================##
  ##   The first phase   ##
  ##---------------------##
  graE<-c()
  graR<-c()
  for(i in seq(1,length(ct),2)){
    if( length(union(graE,graE)) !=length(union(graE,c(ct[i],ct[i+1])))  ){
      graE<-c(graE,ct[i],ct[i+1])
      g<- graph(graE, directed=T)
    }else{
      #print(paste(ct[i],"+",ct[i+1]))
      if(edge_connectivity(g, source = ct[i], target = ct[i+1], checks = TRUE) == 0){
        graE<-c(graE,ct[i],ct[i+1])
        g <- graph(graE, directed=T)
      }
      else{
        graR<-c(graR,ct[i],ct[i+1])
      }
    }
  }

  ##======================##
  ##   The second phase   ##
  ##----------------------##
  print(length(graR))
  g <- graph(graE, directed=F)
  for(i in  seq(1,length(graR),2)){
    #save path as n.x to one_path
    shortpa<-shortest_paths(g, from = graR[i], to = graR[i+1], mode = c("all"))$vpath
    one_path<-names(V(g))[as.integer(shortpa[[1]])]
    #compute cutpoint and MI
    brek <- one_path[-c(1,length(one_path))]
    if(length(brek) == 0){
      #print(" brek=0! ")          ############
      info = 1                     #conscious #
    }else{                         #          #
      #print(" do! ")              ############
      info <- getCMI(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),expMatrix)
      #print(info)
    }
    print(i)  #2016.8.8
    if(info > weight2){ graE <- c(graE,graR[i],graR[i+1]) } #如果互信息大于weight2则添加
  }
  rm("one_path","i","shortpa","info","brek","graR")

  ##=====================##
  ##   The Third phase   ##
  ##---------------------##
  g<- graph(graE, directed=F)
  print(length(graE))
  for(i in  seq(1,length(graE),2)){
    g_d <- g - edge(paste(graE[i],"|",graE[i+1],sep = ""))
    if( edge_connectivity(g_d, source = graE[i], target = graE[i+1], checks = TRUE) > 0){
      shortpa<-shortest_paths(g_d, from = graE[i], to = graE[i+1], mode = c("all"))$vpath
      one_path<-names(V(g))[as.integer(shortpa[[1]])]
      brek <- one_path[-c(1,length(one_path))]
      if(length(brek) > 0){
        print(brek)
        info <- getCMI(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),expMatrix)
        #print(info)
        if(info < weight3){ g <- g_d } #如果互信息小于weight3，那么久确定删除
      }else{
        print(i)
      }
    }
  }
  result <- as_edgelist(g, names = TRUE) # TPDA result
  for(i in 1:length(result[,1])){
    result[i,] <- Translat(result[i,])
  }
  return(data.frame(result))
}



##' Pre-process data to better learn Bayesian networks
##'
##' @param df a data frame containing the the variables in the model
##' @param discrete
##' @return a matrix containing cleand values
##' @importFrom bnlearn discretize
##' @export
dfProcess <- function(df, discrete=F){
  rownames(df) <- df[,1]
  ## remove df`s rowname and colname
  df <- df[,-1]
  ## numerical data frame, for next steps
  tmp <- data.frame(matrix(as.numeric(unlist(df)),ncol = length(df[1,])))
  rownames(tmp) <- rownames(df)
  colnames(tmp) <- colnames(df)
  if(discrete){
  ## dicretize data frame
  ## transfer df
  transdf <- data.frame(t(df))
  tpN <- discretize(transdf,method ='quantile',breaks=2 ) #离散化 三阶段不可离散
  rownames(tpN) <- colnames(df)
  colnames(tpN) <- rownames(df)
  tpN <- data.frame(tpN)
  res <- t(tpN)
  }
  res <- tmp
  return(res)
}

##' Calculating the value of covariance matrix
##'
##' @param vectorA
##' @param vectorB
##' @return
##' @description Calculating the covariance matrix value of vector A and vector B
getDoubleHang_value<-function(vectorA,vectorB){
  Dx1=var(vectorA)
  Dx2=var(vectorB)
  res=Dx1*Dx2-cov(vectorA,vectorB)*cov(vectorA,vectorB)
  return(res)
}

##' Calculating mutual info
##'
##' @param df the file path which you input
##' @param weight
##' @return will be saved as [from,to,value]
##' @export
getMI <- function(df, weight){
  mylist_a<-c()
  mylist_b<-c()
  mylist_info<-c()
  for(xi in 1:nrow(df)){
    mylist_a<-as.double(c(df[xi,]))
    Dxa=var(mylist_a)
    for(yi in 1:nrow(df)){
      #print(yi)
      mylist_b<-as.double(c(df[yi,]))
      Dxb=var(mylist_b)
      cov_ab=getDoubleHang_value(mylist_a,mylist_b)
      mi=0.5*log(Dxa*Dxb/cov_ab)
      mylist_info<-c(mylist_info,mi)
    }
  }
  info<-matrix(mylist_info,nrow(df),nrow(df),byrow = TRUE)

  rname<-c(row.names(df))
  row.names(info)<-row.names(df)
  colnames(info)<-row.names(df)
  # zhao guan xi
  relation<-data.frame()
  for(xi in 1:(nrow(info)-1)){
    yi=xi
    while(TRUE){
      yi=yi+1
      if(nrow(relation)==0){
        from<-c(rname[xi])
        to<-c(rname[yi])
        value<-c(info[xi,yi])
        relation<-data.frame(from,to,value,stringsAsFactors = FALSE)
        relation<-data.frame(xi=c(t(relation)),stringsAsFactors = FALSE)
      }else{
        relation<-data.frame(relation,xi=c(rname[xi],rname[yi],info[xi,yi]),stringsAsFactors = FALSE)
      }
      print(yi)
      if(yi==nrow(info)) break;
    }
  }
  row.names(relation)<-c("from","to","value")
  relation<-data.frame(t(relation),stringsAsFactors = FALSE)

  AllResult<-function(){
    from<-relation[,1][order(as.numeric(relation$value),decreasing = T)]
    to<-relation[,2][order(as.numeric(relation$value),decreasing = T)]
    value<-relation[,3][order(as.numeric(relation$value),decreasing = T)]
    relation2<-data.frame(from,to,value)
    mark1 <- 1
    for(wx in 1:(nrow(relation2))){
      if( as.numeric(as.character(relation2[wx,3])) > weight ){
        print(paste(as.character(relation2[wx,3]),weight))
        mark1 = mark1 +1
      }
    }
    result<-data.frame()
    result<-relation2[1:mark1,]
    return(result)
  }

  tempRst<-AllResult() #函数调用
}


##' Calculating conditional mutual info
##'
##' @param gxi source node
##' @param gyi target node
##' @param cutset cutoff
##' @param gen all genes
##' @return will be saved as [from,to,value]
##' @export
getCMI <- function(gxi, gyi,cutset,gen){
  cutset<-data.frame(cutset,stringsAsFactors = F)
  x<-c()
  y<-c()
  cutset_table<-data.frame()
  for(xi in 1:nrow(gen)){
    if(gxi==gen[xi,1]){
      x<-as.numeric(c(gen[xi,2:length(gen)]))
    }
    if(gyi==gen[xi,]){
      y<-as.numeric(gen[xi,2:length(gen)])
    }
    for(yi in 1:nrow(cutset)){
      if(gen[xi,1]==cutset[yi,1]){
        if(nrow(cutset_table)==0)
          cutset_table<-data.frame(t(gen[xi,]),stringsAsFactors = F)
        else
          cutset_table<-data.frame(cutset_table,t(gen[xi,]),stringsAsFactors = F)
      }
    }
  }
  cutset_table<-data.frame(t(cutset_table),stringsAsFactors = F)
  cutset_table<-cutset_table[,-1]
  to_list<-c()
  cut<-data.frame()
  for(xi in 1:nrow(cutset_table)){
    to_list<-as.numeric(c(cutset_table[xi,]))
    if(nrow(cut)==0)
      cut<-data.frame(to_list)
    else
      cut<-data.frame(cut,to_list)
  }

  if(length(x)==0){
    cut_x<-data.frame(cut)
  }else{
    cut_x<-data.frame(cut,x)
  }
  if(length(y)==0){
    cut_y<-data.frame(cut)
  }else{
    cut_y<-data.frame(cut,y)
  }
  if(length(x)==0 && length(y)==0 ){
    cut_x_y<-data.frame(cut)
  }else if(length(x)==0){
    cut_x_y<-data.frame(cut,y)
  }else if(length(y)==0){
    cut_x_y<-data.frame(cut,x)
  }else{
    cut_x_y<-data.frame(cut,x,y)
  }
  t=det(cov(cut_x))*det(cov(cut_y))/(det(cov(cut))*det(cov(cut_x_y)))
  cmi=0.5*log(t)
}





