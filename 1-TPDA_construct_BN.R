#################################################################
## Using advanced TPDA algorithm to construct Bayesian Network ##
#################################################################
rm(list = ls())

########################
## library R packages ##
########################
library(igraph)
library(bnlearn)

#########################
## setting  thresholds ##
#########################
weight_1 <- 0.01     #First phase MI threshold
weight_2 <- 0.01     #Second phase MI threshold
weight_3 <- 0.01     #Third phase MI threshold


##########################
##  input files` paths  ##
##########################
input_path <- "E:\\GitHub\\TPDA\\testData\\bayesinput.csv"
gene2MI_path <- "E:\\GitHub\\TPDA\\outputData\\1_gene2MI.csv"
network_structure_path <- "E:\\GitHub\\TPDA\\outputData\\2_TDPA_structure.csv"

############################
##  Preprocess Functions  ##
############################
##(1) data preprocessing  
dfProcess <- function(input_path){
  
  zp<-read.csv(input_path,header = F,stringsAsFactors = FALSE)[-1,] ######### 此处输入数据
  row.names(zp)<-zp[,1]
  zp<-zp[-1]
  mydata<-data.frame(zp)  
  tq <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
  rownames(tq) <- rownames(mydata)
  colnames(tq) <- colnames(mydata)
  mydata <- tq
  mydata<-data.frame(t(mydata))
  #tpN <- discretize(mydata,method ='quantile',breaks=2 ) #离散化 三阶段不可离散
  #rownames(tpN) <- rownames(mydata)
  #colnames(tpN) <- colnames(mydata)
  #mydata <- data.frame((tpN))
  z <- t(mydata) ####  z
  return(z)
}

##(2) Calculating the value of covariance matrix 
getDoubleHang_value<-function(mylist1,mylist2){
  Dx1=var(mylist1)
  Dx2=var(mylist2)
  double_h_value=Dx1*Dx2-cov(mylist1,mylist2)*cov(mylist1,mylist2)
}

##(3) Calculating mutual info
## results will be saved as [from,to,value]
getMI <- function(df, weight){
  
  mylist_a<-c()
  mylist_b<-c()
  mylist_info<-c()
  for(xi in 1:nrow(z)){
    mylist_a<-as.double(c(z[xi,]))
    Dxa=var(mylist_a)
    for(yi in 1:nrow(z)){
      #print(yi)
      mylist_b<-as.double(c(z[yi,]))
      Dxb=var(mylist_b)
      cov_ab=getDoubleHang_value(mylist_a,mylist_b)
      mi=0.5*log(Dxa*Dxb/cov_ab)
      mylist_info<-c(mylist_info,mi)
    }
  }
  info<-matrix(mylist_info,nrow(z),nrow(z),byrow = TRUE)
  
  rname<-c(row.names(z))
  row.names(info)<-row.names(z)
  colnames(info)<-row.names(z)
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
      }else
        relation<-data.frame(relation,xi=c(rname[xi],rname[yi],info[xi,yi]),stringsAsFactors = FALSE)
      print(yi)
      if(yi==nrow(info)) 
        break
    }
  }
  
  row.names(relation)<-c("from","to","value")
  relation<-data.frame(t(relation),stringsAsFactors = FALSE)
  
  AllResult<-function(){
    from<-relation[,1][order(as.numeric(relation$value),decreasing = T)] ##降序排序
    to<-relation[,2][order(as.numeric(relation$value),decreasing = T)]
    value<-relation[,3][order(as.numeric(relation$value),decreasing = T)]
    relation2<-data.frame(from,to,value)
    mark1 <- 1
    for(wx in 1:(nrow(relation2))){
      if( as.numeric(as.character(relation2[wx,3])) > weight_1 ){
        print(paste(as.character(relation2[wx,3]),weight_1))
        mark1 = mark1 +1 
      }
    }
    result<-data.frame()
    result<-relation2[1:mark1,]
    return(result)
  } 
  
  tempRst<-AllResult() #函数调用
  
}


##(4) Calculating conditional mutual info
#######################
# prameter description
# gxi source node
# gyi target node
# cutset cutoff
# gen all genes
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


##(5) Three Phase Development Algorithm
#######################
# prameter description
# weight_1,weight_2,weight_3 : mutual info threshold
# input_path : the file path which you input
# gene2MI_path : the file path for outputting gene-gene mutual information
TPDA_algorithm <- function(weight_1,weight_2,weight_3,
                           input_path,
                           gene2MI_path){
  
  gene<-read.csv(input_path,header = FALSE,stringsAsFactors = FALSE)############此处输入数据
  #build_data <- gene2MI ## MI values data table
  build_data <- read.csv(gene2MI_path,header = T,stringsAsFactors = F)
  
  ck <- union(build_data$from,build_data$to)
  lname<-gene[,1]
  gene<-gene[-1,]
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
  snum<-length(build_data$from)
  for(i in 1:snum){ ct<-c(ct,build_data[i,1],build_data[i,2]) }
  for(i in 1:length(ct)){ 
    #print(paste(i,ct[i]))
    ct[i]<-paste("n.",which(labe == ct[i]),sep = "") 
  }
  
  
  ##(4)Start TPDA 
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
      #print(" brek=0! ")          ################      
      info = 1                     #     出错     #
    }else{                         #     位置     #
      #print(" do! ")              ################
      info <- getCMI(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),gene)
      #print(info) 
    }
    print(i)  #2016.8.8
    if(info > weight_2){ graE <- c(graE,graR[i],graR[i+1]) } #如果互信息大于weight_2则添加
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
        info <- getCMI(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),gene)
        #print(info)
        if(info < weight_3){ g <- g_d } #如果互信息小于weight_3，那么久确定删除
      }else{
        print(i)
      }
    }
  }
  result <- as_edgelist(g, names = TRUE) # TPDA result
  for(i in 1:length(result[,1])){
    result[i,] <- Translat(result[i,])
  }
  
}


##########################################
## begin to construct bayesian network ###
##########################################
##(1) data process
z <- dfProcess(input_path)
##(2) get primary mutual information
gene2MI <- getMI(z, weight_1)
write.csv(gene2MI,gene2MI_path,row.names = F)
##(3) start TPDA
result <- TPDA_algorithm(weight_1,weight_2,weight_3,
                         input_path,
                         gene2MI_path,
                         network_structure_path)
##(4) start TPDA
# network_structure_path : the file path for outputting bayesian network structure constructed by TPDA
write.csv(result,network_structure_path,row.names = F)  ## output TPDA result


