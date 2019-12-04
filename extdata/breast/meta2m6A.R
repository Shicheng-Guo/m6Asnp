 
args = commandArgs(trailingOnly=TRUE)

meta.file=as.character(args[1])
output.file=as.character(args[2])

# meta.file="oncoarray_bcac_public_release_oct17.chr22.rsid.txt"
# output.file="oncoarray_bcac_public_release_oct17.chr22.m6A.txt"
 
input<-read.table(meta.file,head=T,sep="\t")
table<-data.frame(input,OR=exp(as.numeric(input$bcac_onco_icogs_gwas_beta)))
db<-read.table("~/hpc/db/hg19/m6ASNP.db147.hg19.uni.txt",head=T)
head(db)
newdb<-db[db$Rs_ID %in% table$phase3_1kg_id,]
out<-data.frame(newdb,table[match(newdb$Rs_ID,table$phase3_1kg_id),])
head(out)
write.table(out,file=output.file,sep="\t",quote=F,row.names =F,col.names = T)
