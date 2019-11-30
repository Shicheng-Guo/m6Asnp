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
