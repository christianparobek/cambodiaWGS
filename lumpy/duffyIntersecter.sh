## 15 December 2015
## To intersect the lumpy output with a bed of DBP location
## I lowered the left limit of DBP coord by 10K to detect the beginning of the dup

for file in `ls sv/`
do

echo $file
echo $file >> lumpy_duffy.txt
bedtools intersect -a sv/$file -b duffy.bed >> lumpy_duffy.txt


done
