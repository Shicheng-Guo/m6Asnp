source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("meta")
library("metafor")
library("survival")
library("survminer")

db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
Symbol2ENSG<-function(Symbol,db){
  ENSG<-as.character(db[match(Symbol,db$V4),8])
  ENSG<-na.omit(data.frame(Symbol,ENSG))
  return(ENSG)
}
ENSG2Symbol<-function(ENSG,db){
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}
ensg2bed<-function(ENSG,db){
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  bed<-unique(db[db$V5 %in% as.character(ENSG),c(1,2,3,5)])
  return(bed)
}

chr2num<-function(x){
x<-output$V1
x<-gsub("chr","",x)
x[x=="X"]<-23
x[x=="Y"]<-24
return(x)
}

chr2num<-function(x){
x<-output$V1
x<-gsub("chr","",x)
x[x=="X"]<-23
x[x=="Y"]<-24
return(x)
}

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")

TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)

# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]

idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
noise<-abs(matrix(rnorm(ncol(input)*nrow(input),0,0.01),nrow=nrow(input),ncol=ncol(input)))
input<-input+noise

input<-log(input+1,2)
input<-RawNARemove(input)
#input<-RawZeroRemove(input)

# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/tsg.positivecontrol.txt",head=F)
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/esophageal/master/phase1.genelist.txt",head=F)
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/breast/master/target.txt",head=T)
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholangiocarcinoma.hg19.bed",head=T,as.is=T,sep="\t")
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/drugtarget/master/extdata/2993drugtarget.txt",sep="\t",as.is=T)[,1]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/encode/master/TFBS/685TFBS.txt",sep="\t",as.is=T)[,1]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cholangiocarcinoma/master/cholangiocarcinoma.V2.hg19.bed",sep="\t",as.is=T)[,4]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cldn6/master/cldn6.kegg.list.txt",sep="\t",as.is=T)[,3]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/lncRNA/master/extdata/tcgameta/1668.ncRNA.dge.txt",sep="\t",as.is=T)[,1]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/encode/master/TFBS/1639TF.txt",sep="\t",as.is=T)[,1]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/epifactors/master/epifactors.txt",sep="\t",as.is=T)[,1]
# xxxv<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/ferroptosis.genelist.csv",head=F)[,2]
# xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/PANC754/coexpression/genelist.txt",sep="\t",as.is=T)[,1]
xxxv<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/m6Asnp/master/m6AeQTLmeta/m6A-eQTL-pick.txt",sep="\t",as.is=T)[,1]

# mkdir /mnt/bigdata/Genetic/Projects/shg047/meta/tfbs/
# mkdir /mnt/bigdata/Genetic/Projects/shg047/meta/tfbs/os
# mkdir /mnt/bigdata/Genetic/Projects/shg047/meta/tfbs/dge

# setwd("~/hpc/methylation/chol/meta/dge")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/drug/dge")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/drug")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/tfbs/os")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/tfbs/dge")
# setwd("/home/guosa/hpc/methylation/chol/meta/dge")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/cldn6")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/lncRNA")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/tfbs")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/epifactor")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/ferroptosis")
# setwd("/mnt/bigdata/Genetic/Projects/shg047/meta/panc754")

# library(limma)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# tab <- getGeneKEGGLinks(species="hsa")
# tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,column="SYMBOL", keytype="ENTREZID")
# head(tab)
# kegg<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/cldn6/master/cldn6.kegg.txt",as.is=T)[,1]
# out<-unique(tab[tab$PathwayID %in% paste("path:",kegg,sep=""),])
# write.table(out,file="cldn6.kegg.list.txt",quote=F,col.names=T,row.names=F,sep="\t")
# write.csv(out,file="cldn6.kegg.list.csv",quote=F)
# xxxv<-as.character(xout$Symbol)

ENSG<-Symbol2ENSG(unique(as.character(xxxv)),db)
xgene<-c(as.character(ENSG[,2]))
ii<-na.omit(unique(unlist(lapply(xgene,function(x) grep(x,rownames(input))))))

Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()
z<-1
for(i in ii){
  z<-z+1
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
  m<-metagen(yi,seTE=vi,data = es,
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=F,
             sm="SMD")

  Symbol<-ENSG2Symbol(rownames(input)[i],db)
  print(c(z,i,as.character(Symbol)))
  pdf(paste(Symbol,"-",rownames(input)[i],".SMD.PANC.pdf",sep=""))
  forest(m,leftlabs = Source,
         lab.e = "Intervention",
         pooled.totals = FALSE,
         smlab = "",studlab=Source,
         text.random = "Overall effect",
         print.tau2 = FALSE,
         col.diamond = "blue",
         col.diamond.lines = "black",
         col.predict = "red",
         print.I2.ci = TRUE,
         digits.sd = 2,fontsize=8)
  dev.off()
}
rownames(rlt)<-rownames(input)[ii]
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt2<-data.frame(rlt)
rlt2<-rlt2[order(rlt2$pval),]
rlt2$symbol<-as.character(ENSG2Symbol(as.character(rownames(rlt2))),db)
head(rlt2)

memo="m6A.dge"

# save dge to csv
write.table(rlt2,file=paste(memo,".tcga.pancancer.smd.meta.pvalue.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)
write.csv(rlt2,file=paste(memo,".tcga.pancancer.smd.meta.pvalue.csv",sep=""),quote=F)

dim(subset(rlt2,pval<10^-8))
up<-(subset(rlt2,beta>0 & pval<10^-8))
down<-(subset(rlt2,beta<0 & pval<10^-8))
write.csv(up,file=paste(memo,".up.tcga.pancancer.smd.meta.pvalue.csv",sep=""),quote=F)
write.csv(down,file=paste(memo,".down.tcga.pancancer.smd.meta.pvalue.csv",sep=""),quote=F)

dir.create("pick")
for(i in 1:80){
x<-paste("cp *-",rownames(rlt2)[i],"*SMD*pdf ./pick",sep="")
system(x)
}

# install.packages("CMplot")
library("CMplot")
ensg<-read.table("~/hpc/db/hg19/ENSG.hg19.bed")
output<-data.frame(ensg[match(unlist(lapply(strsplit(rownames(rlt2),"[.]"),function(x) x[1])),ensg[,5]),],rlt2)

cminput<-na.omit(data.frame(SNP=output$V5,Chromosome=chr2num(output$V1),Position=output$V2,trait1=output$pval,stringsAsFactors=F))
CMplot(cminput,plot.type="b",memo=paste(memo,".fix",sep=""),LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.table(cminput,file=paste(memo,".pval.manhattan.qqplot.meta.dge.1.txt",sep=""),sep="\t",quote=F,row.name=T,col.names=NA)
cminput2<-cminput
cminput2$symbol<-as.character(ENSG2Symbol(as.character(cminput2$SNP),db))
write.table(cminput2,file=paste(memo,".pval.manhattan.qqplot.meta.dge.2.txt",sep=""),sep="\t",quote=F,row.name=T,col.names=NA)
write.csv(cminput2,file=paste(memo,".pval.manhattan.qqplot.meta.dge.2.csv",sep=""),quote=F)



