m6A<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")

brcagwas1<-read.csv("cminput.1.csv",as.is=T)
brcagwas2<-read.csv("cminput.2.csv",as.is=T)
brcagwas3<-read.csv("cminput.3.csv",as.is=T)

out1<-merge(brcagwas1,m6A,by.x="SNP",by.y="V9")
out2<-merge(brcagwas2,m6A,by.x="SNP",by.y="V9")
out3<-merge(brcagwas3,m6A,by.x="SNP",by.y="V9")

rlt1<-subset(out1,as.numeric(trait1)<0.01)
rlt2<-subset(out2,as.numeric(trait1)<0.01)
rlt3<-subset(out3,as.numeric(trait1)<0.01)

rlt<-rbind(rlt1,rlt2,rlt3)

write.csv(rlt,file="brcaGwas.m6A.eQTL.csv",quote=F)
write.table(rlt,file="brcaGwas.m6A.eQTL.txt",quote=F,col.names = NA,row.names = T,sep="\t")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
symbolist<-as.character(unique(rlt[,ncol(rlt)]))
pancancermetadge(symbolist,memo="brca.m6A.eQTL.pancancer.dge")
pancancermetaOsHr(symbolist,memo="brca.m6A.eQTL.pancancer.hr.os")


# esophagus cancer
setwd("/home/guosa/hpc/project/m6A/esophagus")
m6A1<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Esophagus_Gastroesophageal_Junction.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")
m6A2<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Esophagus_Mucosa.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")
m6A3<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Esophagus_Muscularis.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")
rlt<-rbind(m6A1,m6A2,m6A3)

write.csv(rlt,file="esophagus.m6A.eQTL.csv",quote=F)
write.table(rlt,file="esophagus.m6A.eQTL.txt",quote=F,col.names = NA,row.names = T,sep="\t")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
symbolist<-as.character(unique(rlt[,ncol(rlt)]))
pancancermetadge(symbolist,memo="esophagus.m6A.eQTL.pancancer.dge")
pancancermetaOsHr(symbolist,memo="esophagus.m6A.eQTL.pancancer.hr.os")

# liver cancer
wdir<-"/home/guosa/hpc/project/m6A/hcc"
dir.create(wdir)
setwd(wdir)
m6A1<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")
rlt<-rbind(m6A1)
write.csv(rlt,file="esophagus.m6A.eQTL.csv",quote=F)
write.table(rlt,file="esophagus.m6A.eQTL.txt",quote=F,col.names = NA,row.names = T,sep="\t")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
symbolist<-as.character(unique(rlt[,ncol(rlt)]))
pancancermetadge(symbolist,memo="HCC.m6A.eQTL.pancancer.dge")
pancancermetaOsHr(symbolist,memo="HCC.m6A.eQTL.pancancer.hr.os")


# thyroid cancer
setwd("/home/guosa/hpc/project/m6A/")
wdir<-"/home/guosa/hpc/project/m6A/hcc"
dir.create(wdir)
setwd(wdir)
m6A1<-read.csv("/home/guosa/hpc/project/m6A/GTEx_Analysis_v8_eQTL/Thyroid.v8.signif_variant_gene_pairs.txt.rsid.txt.m6A.txt.pick.csv")
rlt<-rbind(m6A1)
write.csv(rlt,file="thyroid.m6A.eQTL.csv",quote=F)
write.table(rlt,file="thyroid.m6A.eQTL.txt",quote=F,col.names = NA,row.names = T,sep="\t")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetadge.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/pancancermetaOsHr.R")
symbolist<-as.character(unique(rlt[,ncol(rlt)]))
pancancermetadge(symbolist,memo="thyroid.m6A.eQTL.pancancer.dge")
pancancermetaOsHr(symbolist,memo="thyroid.m6A.eQTL.pancancer.hr.os")






