###若将三阶段算法结果作为推理的对象，
##那么在读入原始数据时，需要在原始数据中删除三阶段算法中被淘汰的孤立节点【他和其他节点没有边】
rm(list = ls())
setwd("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\Second_qq后的数据\\Gala\\")

################
## 导入程序包 ##
################

library(igraph)
library(bnlearn)
library(Rgraphviz) 


## 设置互信息阈值
weight_1 <- 0.01     #第一阶段互信息阀值
weight_2 <- 0.01     #第二阶段互信息阀值
weight_3 <- 0.01     #第三阶段互信息阀值


## 文件输入路径
pheno_name <- read.csv("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\Second_qq后的数据\\Gala_phenoname.csv", stringsAsFactors = F)
allfiles <- list.files(getwd())

for(kkk in 1:length(allfiles)){
  print(kkk)

input_path <- allfiles[kkk]
c2link_path <- paste0("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\输出数据\\tpda\\"
                              ,strsplit(allfiles[kkk],"qq")[[1]][1],"c2link.csv")
network_structure_path <- paste0("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\输出数据\\tpda\\"
                                         ,strsplit(allfiles[kkk],"qq")[[1]][1],"Structure.csv")
delete_col_path <- paste0("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\输出数据\\tpda\\"
                          ,strsplit(allfiles[kkk],"qq")[[1]][1],"已处理.csv")
output_lisan_path <- paste0(strsplit(allfiles[kkk],"qq")[[1]][1],"_lisan.csv")
pheno_name <- read.csv("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\Second_qq后的数据\\Gala_phenoname.csv", stringsAsFactors = F)

####################
## 数据标准化函数 ## dataframe qqnorm (dfQqnorm)
####################

dfQqnorm <- function(df){

rownames(df) <- df[,1]
df <- df[,-1]
exp2<-df
for(n in 1:length(df[,1])){
  print(n)
  exp2[n,]=as.data.frame(t(as.matrix(qqnorm(df[n,])$x)))
}
exp2 <- cbind(rownames(exp2),exp2)
colnames(exp2)[1] <- "ID"
return(exp2)
}

# e.g 
# df <- read.csv("oil_去掉SNP.csv",stringsAsFactors = F)
# dfpheno <- read.csv("oil_pheno.csv",stringsAsFactors = F)
# dfnormRes <- dfQqnorm(df)
# dfphenormRes <- dfQqnorm(dfpheno)
# write.csv(RFO_qq,"RFObayes_qq.csv",row.names = F)
#df <- read.csv(input_path,stringsAsFactors = F)
#RFO_qq <- dfQqnorm(df)

################################
## 数据格式处理函数dfProcess  ##
################################

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

z <- dfProcess(input_path)
# z <- dfProcess("oilData_qq.csv")

######################################
## 计算列的 协方差矩阵的 行列式的值 ##
######################################

getDoubleHang_value<-function(mylist1,mylist2){
  Dx1=var(mylist1)
  Dx2=var(mylist2)
  double_h_value=Dx1*Dx2-cov(mylist1,mylist2)*cov(mylist1,mylist2)
}



#################
##  互信息计算 ##
#################
## 计算初始互信息，最后的结果relation即三元组[from,to,value]
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

 c2link <- getMI(z, weight_1)
 write.csv(c2link,c2link_path,row.names = F) 



#######################################################
## gene是原始基因，build_data是根据互信息计算的排名  ##
#######################################################

gene<-read.csv(input_path,header = FALSE,stringsAsFactors = FALSE)############此处输入数据
build_data <- read.csv(c2link_path,header = T,stringsAsFactors = F)
ck <- union(build_data$from,build_data$to)
lname<-gene[,1]
gene<-gene[-1,]
labe<-intersect(lname,ck)

#实现将一个n.x转化成GRMZM2G151227
Translat <- function(data_nx){
  retdata <- data_nx
  for(i in 1:length(data_nx)){
    retdata[i] <- labe[which( trans_labe == data_nx[i])]
  }
  return(retdata)
}
#名字转化: GRMZM2G162755类 转为 n.x类
trans_labe<-c()
for(i in 1:length(labe)){ trans_labe<-c(trans_labe,paste("n.",i,sep = "")) }
ct<-c()
snum<-length(build_data$from)
for(i in 1:snum){ ct<-c(ct,build_data[i,1],build_data[i,2]) }
for(i in 1:length(ct)){ 
  #print(paste(i,ct[i]))
  ct[i]<-paste("n.",which(labe == ct[i]),sep = "") 
}


##################
##  三阶段计算  ##
##################

##==============##
##   第一阶段   ##
##--------------##
  
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

#plot(g,layout=layout.fruchterman.reingold, vertex.size=10,vertex.color="green")

##==============##
##   第二阶段   ##
##--------------##

################################################
## 计算条件互信息，源点，目标，割集，所有数据 ##
################################################
Info <- function(gxi,gyi,cutset,gen){
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
  
  #代码翻译计算公式
  #实际计算中有零项，此处区别对待之
  
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

print(length(graR))
g <- graph(graE, directed=F)
for(i in  seq(1,length(graR),2)){
  #找到路径，并存储为n.x到one_path中
  shortpa<-shortest_paths(g, from = graR[i], to = graR[i+1], mode = c("all"))$vpath
  one_path<-names(V(g))[as.integer(shortpa[[1]])]
  #输入n.x输出GRMZM2G162755  for(j in 1:length(one_path)){ one_path[j] <- labe[which( trans_labe == one_path[j])]}
  #计算割点以及互信息
  brek <- one_path[-c(1,length(one_path))]
  if(length(brek) == 0){
    #print(" brek=0! ")          ################      
    info = 1                     #     出错     #
  }else{                         #     位置     #
    #print(" do! ")              ################
    info <- Info(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),gene)
    #print(info) 
  }
  print(i)  #2016.8.8
  if(info > weight_2){ graE <- c(graE,graR[i],graR[i+1]) } #如果互信息大于weight_2则添加
}
rm("one_path","i","shortpa","info","brek","graR")

##==============##
##   第三阶段   ##
##--------------##
g<- graph(graE, directed=F)
print(length(graE))
#起点和终点，其中每一个相邻的两个节点都代表着一条边，除的关系是(i+1)/2,其中所有的边都没有重复的边
for(i in  seq(1,length(graE),2)){
  g_d <- g - edge(paste(graE[i],"|",graE[i+1],sep = ""))  #删除当前边并存放在g_d中
  #在 g_d 中查找是不是还存在路径,存在的话就进行计算互信息，不存在则跳过
  if( edge_connectivity(g_d, source = graE[i], target = graE[i+1], checks = TRUE) > 0){
    #找到路径，并存储为n.x到one_path中
    shortpa<-shortest_paths(g_d, from = graE[i], to = graE[i+1], mode = c("all"))$vpath
    one_path<-names(V(g))[as.integer(shortpa[[1]])]
    #计算割点以及互信息
    brek <- one_path[-c(1,length(one_path))]
    if(length(brek) > 0){
      print(brek)
      info <- Info(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),gene)
      #print(info)
      if(info < weight_3){ g <- g_d } #如果互信息小于weight_3，那么久确定删除
    }else{
      print(i)
    }
  }
}
#最终结果在result中
result <- as_edgelist(g, names = TRUE)
for(i in 1:length(result[,1])){
  result[i,] <- Translat(result[i,])
}

#result <- threePhase()

write.csv(result,network_structure_path,row.names = F)  ##输出三阶段算法后的结果


##################################
## 对三阶段结果中的节点匹配数据 ##
##################################

infer_gene<-data.frame()
the_origin_input<-read.csv(input_path,header = F)[-1,]
gene_name<-union(result[,1],result[,2])
#infer_gene <- the_origin_input[gene_name %in% the_origin_input[,1],]

for(i in 1:length(gene_name))
{
  flag=0
  for(j in 1:nrow(the_origin_input))
  {
    if(as.character(gene_name[i])==as.character(the_origin_input[j,1]))
    {
      infer_gene<-rbind(infer_gene,the_origin_input[j,])
      flag=1
    }
    if(flag==1)
    {break}
  }
}

write.csv(infer_gene,delete_col_path,row.names=F)

##############################################
##人工处理掉Infer_gene的第一列
#然后供程序读入

input_path_infer_gene <- delete_col_path
zp<-read.csv(input_path_infer_gene,header = F,stringsAsFactors = FALSE)[-1,]
row.names(zp)<-zp[,1]
zp<-zp[-1]
mydata<-data.frame(zp)  
tq <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
rownames(tq) <- rownames(mydata)
colnames(tq) <- colnames(mydata)
mydata <- tq
mydata<-data.frame(t(mydata))

tpN <- discretize(mydata,method ='quantile',breaks=2 ) #离散化
rownames(tpN) <- rownames(mydata)
colnames(tpN) <- colnames(mydata)
mydata <- data.frame((tpN))  


###########################################
##  用三阶段结果和graphNEL函数构建网络   ##
###########################################

library(graph)
result<-read.csv(network_structure_path,stringsAsFactors = FALSE) ################################>>>>将三阶段得到的结果作为输入【两列的csv文件】

v<-union(result[,1],result[,2])
edL=vector("list",length=length(v)) ##这里length(v)为v中所有节点数目
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

gRbn<-as.bn(gR)  ##将gR化为bn格式，才能进行推理



##################对网络进行推理学习，每1个基因计算他对性状的影响概率
###首先得到每个基因的k个不重复离散区间, ex : 2值离散就得到两个区间
###因为三阶段结果是八百多个基因，所以进行离散的也应该是这八百多个基因，


lisan<-data.frame()## 2
for(i in 1:ncol(mydata))
{
  xx<-unique(mydata[,i])
  lisan<-rbind(lisan,t(as.data.frame(xx)))
}


######################################################################文件格式化操作
write.csv(lisan,output_lisan_path)              ##将得到的不重复离散区间写入文件
add_rowname<-read.csv(input_path_infer_gene) 
all_rowname<-add_rowname[,1]                    ##取出基因名称
no_repeat<-read.csv(output_lisan_path)   
no_repeat$X<-add_rowname[,1]                    ##将基因名称和得到的区间整合
write.csv(no_repeat,output_lisan_path)          ##可供下次读入的最终结果

#############################################################################
re<-read.csv(output_lisan_path)
re<-re[1:nrow(re),2:ncol(re)]
############################################################################
gene_1<-c(1)
gene_2<-c(2)



###########################################################################
##infer_result<-data.frame(gene_1,gene_2,type_LUT,type_ZEA,type_AC,type_BC)

for(kk in 1:nrow(pheno_name))
{
  print(kk)
  Marker_name <- as.character(pheno_name[kk,1]) 
  ###########################################################################
  ##infer_result<-data.frame(gene_1,gene_2,type_LUT,type_ZEA,type_AC,type_BC)
  #print_all<-data.frame(ID=c(1:26)) ##此数据框 将每次两两推理得到的8个值 都拿出来
  five_count<-data.frame(ID=c(1:5)) ###  五值的相加
  cursor<-data.frame(ID=c(1:5))  ##满足超过阈值的基因 将被加入此cursor数据框
  
  set.seed(1000)##设置种子点
  fitted_time<-system.time(fitted<-bn.fit(gRbn,data = mydata,method='bayes')) ###这里可以改变，选择是对三阶段的结果推理或者对九种算法中的某种推理
  fitted_time
  ####################################################################################################
  for(i in 1:nrow(re))
  {
    
    if(as.character(re[i,1])==Marker_name)
    { print(i)
      E_str_1=as.character(re[i,2])  
      E_str_2=as.character(re[i,3])  
    }
    
    
  }
  
  
  E_str_type_1=paste("(",Marker_name,"=='",E_str_1,"')" ,sep="")
  E_str_type_2=paste("(",Marker_name,"=='",E_str_2,"')" ,sep="")
  
  
  #############################################################################################
  infer_time<-system.time(
    
    for(i in 1:nrow(re))
    { 
      print_all_query<-c()
      five_count_query<-c()
      node_name<-as.character(re[i,1])  ##第i列的 节点的名称取出来
      pp1<-as.character(re[i,2])
      pp2<-as.character(re[i,3])
      str1=paste(node_name,"=='",pp1,"'",sep="")
      str2=paste(node_name,"=='",pp2,"'",sep="")
      
      x1=cpquery(fitted,event = eval(parse(text=E_str_type_1)),evidence=eval(parse(text=str1)),method="ls")
      x2=cpquery(fitted,event = eval(parse(text=E_str_type_1)),evidence=eval(parse(text=str2)),method="ls")
      x3=cpquery(fitted,event = eval(parse(text=E_str_type_2)),evidence=eval(parse(text=str1)),method="ls")
      x4=cpquery(fitted,event = eval(parse(text=E_str_type_2)),evidence=eval(parse(text=str2)),method="ls")
      
      # print_all_query<-c(node_name,"E3_n1569",as.character(x1),as.character(x2),as.character(x3),as.character(x4))
      
      
      y1=x1
      y2=x2
      y3=x3
      y4=x4
      
      five_count_query<-c(node_name,as.character(y1),as.character(y2),as.character(y3),as.character(y4))
      five_count<-data.frame(five_count,five_count_query)
      
      ##设置阈值，将概率大于某个数的点输出
      # if(y1>=2.55|y2>=2.55|y3>=2.55|y4>=2.55)
      #{
      #  cursor_name<-c(node_name,as.character(y1),as.character(y2),as.character(y3),as.character(y4))
      #  cursor<-data.frame(cursor,cursor_name)  
      # }
      
      ##str11=paste("(",names(mydata)[3],"=='",p1,"')",sep="") 
      ##cpquery(fitted,event=eval(parse(text = str11)),evidence = (LUT=='[-2.93,0]') ,method = "ls")
      #cpquery(fitted,event = (LUT=='(0,2.93]'),evidence =((GRMZM2G304575=='[-2.65,-0.0648]')&(GRMZM2G304573=='[-3,-0.017]')&(GRMZM2G106479
      ##=='[-3,-0.0307]')) ,method = "ls")
      ##lung.sample<-cpdist(fitted,nodes="LUT",evidence =(GRMZM2G012966=='[-3,-0.0511]'),method="ls")
      
    }
  )
  
  #write.csv(print_all,"D://标准数据基因一推一_三阶段.csv")
  output<-paste("E:\\学习资料存放处\\贝叶斯网络\\陈圆圆师姐\\RFO金敏整理给刘老师\\输出数据\\tpda\\"
                ,strsplit(allfiles[kkk],"qq")[[1]][1],Marker_name,"_infer_tpda.csv",sep="")
  
  
  write.csv(five_count,output,row.names = F)
  
}
}
#write.csv(cursor,"E:\\beysian_network\\三阶段_Bcry_811_输出超过阈值的点_2.csv")
