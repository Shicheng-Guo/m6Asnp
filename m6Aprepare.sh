## Download all the m6A-SNP files

## merge all the m6A-SNPs
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

input=m6ASNP.db147.hg19.uni.txt
perl -p -i -e 's/Functional //g' $input
perl -p -i -e 's/' $input
perl -p -i -e 's/miCLIP://g' $input
perl -p -i -e 's/miCLIP\&//g' $input
perl -p -i -e 's/PA-m6A-Seq://g' $input
perl -p -i -e 's/\(//g' $input
perl -p -i -e 's/\)//g' $input
perl -p -i -e 's/ SNV//g' $input
perl -p -i -e 's/MeRIP-Seq://g' $input
perl -p -i -e 's/Prediction://g' $input


## prepare the txt files for m6A-SNP
data<-read.table("m6ASNP.db147.hg19.txt",head=T,sep="\t")
data$start=data$Position-1
data<-data[,c(1,11,2:10)]
data<-unique(data)
write.table(data,file="m6ASNP.db147.hg19.uni.bed",sep="\t",quote=F,col.names = T,row.names = F)

input=m6ASNP.db147.hg19.uni.bed
perl -p -i -e 's/Functional //g' $input
perl -p -i -e 's/' $input
perl -p -i -e 's/miCLIP://g' $input
perl -p -i -e 's/miCLIP\&//g' $input
perl -p -i -e 's/PA-m6A-Seq://g' $input
perl -p -i -e 's/\(//g' $input
perl -p -i -e 's/\)//g' $input
perl -p -i -e 's/ SNV//g' $input
perl -p -i -e 's/MeRIP-Seq://g' $input
perl -p -i -e 's/Prediction://g' $input

