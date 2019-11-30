for i in {1..23} 
do
awk -v pat=$i '$4==pat' meta_v3_onco_euro_overall_ChrAll_1_release.txt > ./chr/meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.txt &
echo $i
done

