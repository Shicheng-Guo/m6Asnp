### Genome-wide identificaiton of m6A-SNP and human disease


```
cd ~/hpc/project/m6A
```
Timeline: 

* 2020/01/10: download GTEx eqtl data to `/home/guosa/hpc/db/dbSNP153/hg38/dbSNP153.GRCh38p12b.vcf.gz`
* 2020/01/10: Pan-cancer GWAS and m6SNP collection and check TCGA experssion database to overlap with DGE
* 2019/10/28: [GWAS-Catlog](https://www.ebi.ac.uk/gwas/docs/file-downloads): 11/27/2019, [161,526 records](rsid.txt), 103,202 unique SNPs. 
* 2018/12/12: m6AVar! Database of functional variants involved in m6A modification: http://m6avar.renlab.org/
* extract all m6A-SNP: `perl -lane '{print $1 if (/(rs\d+)/)}' Human_*.txt | sort -u > m6Asnp.txt`
