### Genome-wide identificaiton of m6A-SNP and human disease

* [GWAS-Catlog](https://www.ebi.ac.uk/gwas/docs/file-downloads): 11/27/2019
* awk -F"\t" '{print $22,$25,$8}' OFS="\t" gwas_catalog_v1.0-associations_e96_r2019-11-21.tsv > rsid.txt

