setwd("/gpfs/home/guosa/hpc/project/m6A")

data<-read.table("gwas_catalog_v1.0-associations_e96_r2019-11-21.tsv",head=T,sep="\t")
data<-read.table("rsid.tsv",head=T,sep="\t")

## input all m6A-SNPs
setwd("/gpfs/home/guosa/hpc/project/m6A")
file=list.files(pattern="Human")
snp<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  rlt<-temp[temp$Rs_ID %in% data$SNPS,]
  print(file[i])
  if(nrow(rlt)<10){
  print(rlt)
  }
  snp<-c(snp,as.character(rlt$Rs_ID))
}
snp<-unique(snp)
input1<-unique(data[data$SNPS %in% snp,])
write.table(input,file="m6A-GWAS-hit-overall.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)

i=4
temp<-read.table(file[i],head=T,sep="\t")
input2<-unique(temp[temp$Rs_ID %in% snp,])

write.table(input2,file="m6A-GWAS-hit-m6AVAR.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)

snp1<-as.character(input1$SNPS)
snp2<-as.character(input2$Rs_ID)

input3<-data.frame(input1,input2[match(snp1,snp2),])
write.table(input3,file="m6A-GWAS-hit-m6AVAR-final.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)


setwd("/gpfs/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL")

SNP<-c()
file=list.files(pattern="*.txt")
for(i in 1:length(file)){
  temp<-read.table(file[i],sep="\t",head=T)
  temp<-subset(temp,qval<0.05)
  SNP<-c(SNP,as.character(temp$rs_id_dbSNP151_GRCh38p7))
}
SNP<-unique(SNP)

input4<-input3[input3$SNPS %in% SNP,]
write.table(input4,file="../m6A-GWAS-hit-m6AVAR-eQTL.20191127.txt",sep="\t",quote=F,col.names=NA,row.names = T)


library("CMplot")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A")
P1<-read.table("m6A-GWAS-hit-alzheimer.stage1.20191127.txt",head=T)
P2<-read.table("m6A-GWAS-hit-alzheimer.stage2.20191127.txt",head=T)
P<-rbind(P1,P2)
head(P)
cminput<-data.frame(SNP=P$MarkerName,Chromosome=P$Chromosome,Position=P$Position,trait1=P$Pvalue)
CMplot(cminput,plot.type="b",ylim=25,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A/prostate")
P<-read.table("meta_v3_onco_euro_overall_ChrAll_1_release.Guo.txt",head=T)
head(P)
cminput<-data.frame(SNP=P$MarkerName,Chromosome=P$Chromosome,Position=P$Position,trait1=P$Pvalue)
CMplot(cminput,plot.type="b",ylim=25,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


table<-subset(P,Pvalue<10^-3)
table<-data.frame(table,OR=exp(table$Beta))
write.table(table,file="Table1.txt",sep="\t",quote=F,row.names =F,col.names = T)


###########################################
###########################################
###########################################
data<-read.table("/gpfs/home/guosa/hpc/project/m6A/prostate/meta_v3_onco_euro_overall_ChrAll_1_release.Guo.txt",head=T,sep="\t")
setwd("/gpfs/home/guosa/hpc/project/m6A")
file=list.files(pattern="Human")
snp<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  rlt<-temp[temp$Rs_ID %in% data$SNP,]
  print(file[i])
  if(nrow(rlt)<10){
    print(rlt)
  }
  snp<-c(snp,as.character(rlt$Rs_ID))
}
snp<-unique(snp)
input<-unique(data[data$SNP %in% snp,])
write.table(input,file="meta_v3_onco_euro_overall_ChrAll_1_release.m6A.txt",sep="\t",quote=F,col.names=NA,row.names = T)

library("CMplot")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A/prostate/chr")

input<-read.table("meta_v3_onco_euro_overall_ChrAll_1_release.m6A.txt",head=T,sep="\t")
cminput<-data.frame(SNP=input$SNP,Chromosome=input$Chr,Position=input$position,trait1=input$Pvalue)
CMplot(cminput,plot.type="b",ylim=25,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


input<-read.table("meta_v3_onco_euro_overall_ChrAll_1_release.m6A.txt",head=T,sep="\t")
table<-subset(input,Pvalue<10^-2)
head(table)
table<-data.frame(table,OR=exp(table$Effect))
head(table)
db<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/m6ASNP.db147.hg19.uni.txt",head=T)
head(db)
newdb<-db[db$Rs_ID %in% table$SNP,]
out<-data.frame(newdb,table[match(newdb$Rs_ID,table$SNP),])
head(out)
write.table(out,file="Table1.txt",sep="\t",quote=F,row.names =F,col.names = T)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A/prostate/chr")
file=list.files(pattern="*.m6A.txt")
input<-c()
for(i in 1:length(file)){
temp<-read.table(file[i],head=T,sep="\t")
input<-rbind(input,temp)
}
head(input)
cminput<-data.frame(SNP=input$SNP,Chromosome=input$Chr,Position=input$position,trait1=input$Pvalue)
head(cminput)
CMplot(cminput,plot.type="b",ylim=25,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

table<-subset(input,Pvalue<10^-2)
head(table)
table<-data.frame(table,OR=exp(table$Effect))
eqtl<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/eQTL/GTEx/GTEx_Analysis_v8_eQTL/eqtl.v8.txt",head=T,sep="\t",row.names = 1)
eqtl[1:5,1:4]
output<-data.frame(table,eqtl[match(table$SNP,rownames(eqtl)),])
write.table(output,file="Table1.txt",sep="\t",quote=F,row.names =F,col.names = T)
######################################################################################
setwd("/gpfs/home/guosa/hpc/project/m6A")
file=list.files(pattern="Human")
tR<-c("Chromosome","Position","Rs_ID","Reference_Base","Alterative_Base","Gene","MutType","m6A_Function","MutRegion","Confidence_Level")
merge<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  temp<-temp[,match(tR,colnames(temp))]
  merge<-rbind(merge,temp)
  print(i)
}
write.table(merge,file="m6ASNP.db147.hg19.txt",sep="\t",quote=F,col.names = T,row.names = F)


## prepare the bed files for m6A-SNP
data<-read.table("m6ASNP.db147.hg19.txt",head=T,sep="\t")
data$start=data$Position-1
data<-data[,c(1,11,2:10)]
data<-unique(data)
write.table(data,file="m6ASNP.db147.hg19.uni.txt",sep="\t",quote=F,col.names = T,row.names = F)

perl -p -i -e 's/Functional //g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/miCLIP://g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/miCLIP://g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/miCLIP\&//g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/PA-m6A-Seq://g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/\(//g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/\)//g' m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/ SNV//g' m6ASNP.db147.hg19.uni.txt

## prepare the txt files for m6A-SNP
data<-read.table("m6ASNP.db147.hg19.txt",head=T,sep="\t")
data$start=data$Position-1
data<-data[,c(1,11,2:10)]
data<-unique(data)
write.table(data,file="m6ASNP.db147.hg19.uni.bed",sep="\t",quote=F,col.names = T,row.names = F)

perl -p -i -e 's/Functional //g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/miCLIP://g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/miCLIP://g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/miCLIP\&//g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/PA-m6A-Seq://g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/\(//g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/\)//g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/ SNV//g' m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/MeRIP-Seq://g' m6ASNP.db147.hg19.uni.bed





