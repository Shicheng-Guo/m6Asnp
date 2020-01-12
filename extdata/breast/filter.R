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

