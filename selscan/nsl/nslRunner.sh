## A script that will run selscan
## Started 04 December 2015
## Specify the starting VCF and a folder for work
## Usage: bash selscanRunner.sh </path/to/in.vcf> <folder_for_work> <species_pf_or_pv> <svg_pic_name>


#### READ COMMAND LINE ####

vcf=$1
wkdir=$2
species=$3


#### MAKE VCF DIPLOID STYLE ####
sed 's/\t0:/\t0\/0:/g' $vcf | sed 's/\t1:/\t1\/1:/g' > selscan/nsl/diploid.tmp
cd selscan/nsl
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


#### NORMALIZE nSL ####
cd $wkdir
norm --ihs --files chr01.res.nsl.out chr02.res.nsl.out chr03.res.nsl.out chr04.res.nsl.out chr05.res.nsl.out chr06.res.nsl.out chr07.res.nsl.out chr08.res.nsl.out chr09.res.nsl.out chr10.res.nsl.out chr11.res.nsl.out chr12.res.nsl.out chr13.res.nsl.out chr14.res.nsl.out

#### PLOT nSL ####
cd ..
Rscript nsl.r $wkdir $species


