######$=dx

cd /gpfs/home/guosa/hpc/project/m6A/breast

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
echo Rscript meta2m6A.R meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.rsid.txt meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.m6A.txt >> $i.job
qsub $i.job
done
