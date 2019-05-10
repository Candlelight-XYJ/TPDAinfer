
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

gRbn<-as.bn(gR)  ##convert gR to bnlearn object



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


## parameter learning
set.seed(1000)##设置种子点
fitted_time<-system.time(fitted<-bn.fit(gRbn,data = mydata,method='bayes')) ###这里可以改变，选择是对三阶段的结果推理或者对九种算法中的某种推理
fitted_time