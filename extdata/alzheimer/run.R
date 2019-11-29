
# Stage I
setwd("/gpfs/home/guosa/hpc/project/m6A/alzheimer")
data<-read.table("/gpfs/home/guosa/hpc/project/m6A/alzheimer/Kunkle_etal_Stage1_results.txt",head=T)
setwd("/gpfs/home/guosa/hpc/project/m6A")
file=list.files(pattern="Human")
snp<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  rlt<-temp[temp$Rs_ID %in% data$MarkerName,]
  print(file[i])
  if(nrow(rlt)<10){
  print(rlt)
  }
  snp<-c(snp,as.character(rlt$Rs_ID))
}
snp<-unique(snp)
input1<-unique(data[data$MarkerName %in% snp,])
write.table(input1,file="m6A-GWAS-hit-alzheimer.stage1.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)

# Stage II
setwd("/gpfs/home/guosa/hpc/project/m6A/alzheimer")
data<-read.table("/gpfs/home/guosa/hpc/project/m6A/alzheimer/Kunkle_etal_Stage2_results.txt",head=T)
setwd("/gpfs/home/guosa/hpc/project/m6A")
file=list.files(pattern="Human")
snp<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  rlt<-temp[temp$Rs_ID %in% data$MarkerName,]
  print(file[i])
  if(nrow(rlt)<10){
  print(rlt)
  }
  snp<-c(snp,as.character(rlt$Rs_ID))
}
snp<-unique(snp)
input2<-unique(data[data$MarkerName %in% snp,])
write.table(input2,file="m6A-GWAS-hit-alzheimer.stage2.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)

## manhattan plot and qq-plot
library("CMplot")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A")
P1<-read.table("m6A-GWAS-hit-alzheimer.stage1.20191127.txt",head=T)
P2<-read.table("m6A-GWAS-hit-alzheimer.stage2.20191127.txt",head=T)
P<-rbind(P1,P2)
head(P)
cminput<-data.frame(SNP=P$MarkerName,Chromosome=P$Chromosome,Position=P$Position,trait1=P$Pvalue)
CMplot(cminput,plot.type="b",ylim=20,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)



