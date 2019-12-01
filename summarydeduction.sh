
for i in {1..23}
do
awk '{print $1,$3,$12,$13,$14}' meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.rsid.txt | grep rs | grep -v AA | grep -v TT | grep -v TA | grep -v GT | grep -v CA |grep -v AG| grep -v AT| grep -v AC | grep -v TG | grep -v GA > meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.pc.rsid.Guo.txt 
gzip meta_v3_onco_euro_overall_ChrAll_1_release.chr$i.pc.rsid.Guo.txt 
done
