args = commandArgs(trailingOnly=TRUE)

meta.file=as.character(args[1])
output.file=as.character(args[2])

input<-read.table(meta.file,head=T,sep="\t")
table<-subset(input,Pvalue<=1)
table<-data.frame(table,OR=exp(table$Effect))
db<-read.table("~/hpc/db/hg19/m6ASNP.db147.hg19.uni.txt",head=T)
head(db)
newdb<-db[db$Rs_ID %in% table$SNP,]
out<-data.frame(newdb,table[match(newdb$Rs_ID,table$SNP),])
head(out)
write.table(out,file=output.file,sep="\t",quote=F,row.names =F,col.names = T)
