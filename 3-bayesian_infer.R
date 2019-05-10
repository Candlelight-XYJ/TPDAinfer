##################对网络进行推理学习，每1个基因计算他对性状的影响概率
###首先得到每个基因的k个不重复离散区间, ex : 2值离散就得到两个区间

lisan<-data.frame()
for(i in 1:ncol(mydata))
{
  xx<-unique(mydata[,i])
  lisan<-rbind(lisan,t(as.data.frame(xx)))
}

######################################################################文件格式化操作
##write.csv(lisan,output_lisan_path)              ##将得到的不重复离散区间写入文件
add_rowname<-read.csv(input_path) 
all_rowname<-add_rowname[,1]                    ##取出基因名称
##no_repeat<-read.csv(output_lisan_path)   
no_repeat<-lisan
no_repeat$X<-add_rowname[,1]                    ##将基因名称和得到的区间整合
write.csv(no_repeat,output_lisan_path)          ##可供下次读入的最终结果

#############################################################################
##re<-read.csv(output_lisan_path)
re<-no_repeat
gene_1<-c(1)
gene_2<-c(2)
