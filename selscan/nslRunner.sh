## A script that will run selscan
## Started 04 December 2015
## Specify the starting VCF and a folder for work
## Usage: bash selscanRunner.sh </path/to/in.vcf> <folder_for_work> <species_pf_or_pv> <svg_pic_name>


#### READ COMMAND LINE ####

vcf=$1
wkdir=$2
species=$3
picname=$4


#### MAKE VCF DIPLOID STYLE ####
sed 's/\t0:/\t0\/0:/g' $vcf | sed 's/\t1:/\t1\/1:/g' > selscan/diploid.tmp
cd selscan
mkdir -p $wkdir

#### SPLIT VCF BY CHR ####

if [ $species = "pv" ]
then

	for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	do

	grep "^#\|^Pv_Sal1_chr$i" diploid.tmp > $wkdir/chr$i.vcf

	done

elif [ $species = "pf" ]
then

	for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	do

	grep "^#\|^Pf3D7_$i\_v3" diploid.tmp > $wkdir/chr$i.vcf

	done

else
	echo "Something is wrong with your species notation. Choose either pf or pv."

fi

rm diploid.tmp ## cleanup


#### RUN nSL ####

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

selscan --nsl --max-extend-nsl 50 --vcf $wkdir/chr$i.vcf --out $wkdir/chr$i.res

done


#### PLOT nSL ####

Rscript nsl.r $wkdir $species $picname
