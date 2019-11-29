awk '$8<0.00001{print}' Kunkle_etal_Stage1_results.txt > Kunkle_etal.txt
awk '$8<0.00001{print}' Kunkle_etal_Stage2_results.txt >> Kunkle_etal.txt


