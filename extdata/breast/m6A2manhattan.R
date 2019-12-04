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
