for i in {1..23} 
do
head -n 1   oncoarray_bcac_public_release_oct17.txt > ./chr/oncoarray_bcac_public_release_oct17.chr$i.txt
echo $i
done

for i in {1..23} 
do
awk -v pat=$i '$3==pat' oncoarray_bcac_public_release_oct17.txt >> ./chr/oncoarray_bcac_public_release_oct17.chr$i.txt &
echo $i
done

mkdir temp
for i in {1..23}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo perl addrs.pl $i \> oncoarray_bcac_public_release_oct17.chr$i.rsid.txt >> $i.job
qsub $i.job
done



mkdir temp
for i in {1..23}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo Rscript meta2m6A.R oncoarray_bcac_public_release_oct17.chr$i.rsid.txt oncoarray_bcac_public_release_oct17.chr$i.m6A.txt>> $i.job
qsub $i.job
done


# meta2m6A.R
args = commandArgs(trailingOnly=TRUE)
meta.file=as.character(args[1])
output.file=as.character(args[2])
input<-read.table(meta.file,head=T,sep="\t")
table<-data.frame(table,OR=exp(table$bcac_onco_icogs_gwas_beta))
db<-read.table("~/hpc/db/hg19/m6ASNP.db147.hg19.uni.txt",head=T)
head(db)
newdb<-db[db$Rs_ID %in% table$phase3_1kg_id,]
out<-data.frame(newdb,table[match(newdb$Rs_ID,table$phase3_1kg_id),])
head(out)
write.table(out,file=output.file,sep="\t",quote=F,row.names =F,col.names = T)

# m6A2manhattan
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/m6A/breast/chr")
file<-list.files(pattern="*.m6A.txt")
input<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  input<-rbind(input,temp)
}
head(input)
cminput<-data.frame(SNP=input$Rs_ID,Chromosome=input$chr,Position=input$Position,trait1=input$bcac_onco_icogs_gwas_P1df)

CMplot(cminput,plot.type="b",ylim=ymax,LOG10=TRUE,file="jpg",memo="bcac_onco_icogs_gwas_P1df",
       threshold=c(1e-6,1e-4),threshold.lty=c(1,2),signal.col=c("red","green"),dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)

cminput<-data.frame(SNP=input$Rs_ID,Chromosome=input$chr,Position=input$Position,trait1=input$bcac_onco_icogs_gwas_erpos_P1df)
CMplot(cminput,plot.type="b",ylim=ymax,LOG10=TRUE,file="jpg",memo="bcac_onco_icogs_gwas_erpos_P1df",
       threshold=c(1e-6,1e-4),threshold.lty=c(1,2),signal.col=c("red","green"),dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)

cminput<-data.frame(SNP=input$Rs_ID,Chromosome=input$chr,Position=input$Position,trait1=input$bcac_onco_icogs_gwas_erneg_P1df)
CMplot(cminput,plot.type="b",ylim=ymax,LOG10=TRUE,file="jpg",memo="bcac_onco_icogs_gwas_erneg_P1df",
       threshold=c(1e-6,1e-4),threshold.lty=c(1,2),signal.col=c("red","green"),dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6)
