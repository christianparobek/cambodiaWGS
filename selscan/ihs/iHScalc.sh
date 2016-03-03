## A bash script
## Does a number of things:
##	1) Makes likelihood files that are for the right number of individuals
##	2) Calculates genetic distance using LDhat
##	3) Converts those genetic distance files to MAP format files
##	4) Uses the MAP files and VCF files to calculate iHS using selscan
## 
## USAGE: bash iHScalc.sh <input_vcf> <wkdir_for_results>  <species:_pv_or_pf>
##
## example: http://classic.sciencemag.org/content/suppl/2012/03/14/science.1216872.DC1/Auton.SOM.pdf

#############################################
############ PARSE COMMAND LINE #############
#############################################

vcf=$1
dir=$2
species=$3
ehhcutoff=$4

cd selscan/ihs


#############################################
####### CHOOSE CHR FORMAT BY SPECIES ########
#############################################

chrtype=''

if [ $species = "pv" ]
then
	chrtype="Pv_Sal1_chr\$i"
elif [ $species = "pf" ]
then
	chrtype="Pf3D7_\$i\_v3"
else
	echo "Something is wrong with your species notation. Choose either pf or pv."
fi

#############################################
####### FOR LOOP FOR EACH CHROMOSOME ########
#############################################

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do 

chr="$(eval echo $chrtype)"

mkdir -p $dir
	# make directory, if it doesn't already exist

vcftools --vcf $vcf --out $dir/$i --ldhat --phased --chr $chr
	# convert vcf to ldhat format using vcftools
	# outputs in diploid format, so need to fix

dipnum=`head -1 $dir/$i.ldhat.sites | cut -f1` # get number of diploid chromosomes
hapnum=$(($dipnum / 2)) # divide that by two to get number of haploid chrs (i.e. = our num indiv)

grep -v '???????????????\|\-1' $dir/$i.ldhat.sites | sed 's/'"$dipnum\t"'/'"$hapnum\t"'/' | sed 's/-0//' > $dir/tmp_$i
	# convert ldhat sites file from diploid format to haploid
	# remove the extra chrs with question marks
	# halve the nubmer of individuals listed in the header
rm $dir/$i.ldhat.sites # remove the diploid ldhat.sites file
mv $dir/tmp_$i $dir/$i.ldhat.sites # rename the haploid ldhat.sites file


	######### RUN LDhat lkgen, IF NEEDED #########

if ! [ -e lk_files/lk_n$hapnum\_t0.001 ] # test to see if correct size lookup file exists
then
	lkgen -lk lk_files/lk_n100_t0.001 -nseq $hapnum
		# going to use Watterson's theta - 0.001 because that's roughly what Chang 2012 estimated
	mv new_lk.txt lk_files/lk_n$hapnum\_t0.001
fi

	######### RUN LDhat interval #########

interval -seq $dir/$i.ldhat.sites -loc $dir/$i.ldhat.locs -lk lk_files/lk_n$hapnum\_t0.001 -its 100000 -bpen 5 -samp 2000
	# -bpen should be optimized in range 0-50 to balance sens and spec; maybe use MCMC
	# want to use 10000000 its, but testing a smaller number
mv rates.txt $dir/$i.rates # rename the output rates file

	######### RUN LDhat stat #########

/proj/julianog/bin/LDhat-2.2/stat -input $dir/$i.rates -burn 5 -loc $dir/$i.ldhat.locs
	# look at stat to see what's cookin'
	# I think this excludes the first 1,000,000 iterations as burn in?
	# Manual recommends to exclude the first 100K as burn in
	# was set to 500 but Im setting to 10 to accommodate the fewer its
mv res.txt $dir/$i.res # rename the results file

# can maybe assess convergence using the provided R scripts (manual p23)
# Do I need to assess covnergence over a range of block penalties?

	######### CONVERT res TO MAP #########
Rscript res2map.r --file=$dir/$i.res --chr=$chr --out=$dir/$i.map

	######### CONVERT res TO MAP #########
selscan --ihs --vcf /proj/julianog/users/ChristianP/cambodiaWGS/selscan/$dir/chr$i.vcf --map $dir/$i.map --out $dir/$i --cutoff $ehhcutoff

done

#############################################
####### NORMALIZE ACROSS CHROMOSOMES ########
#############################################

cd $dir
norm --ihs --files 01.ihs.out 02.ihs.out 03.ihs.out 04.ihs.out 05.ihs.out 06.ihs.out 07.ihs.out 08.ihs.out 09.ihs.out 10.ihs.out 11.ihs.out 12.ihs.out 13.ihs.out 14.ihs.out

cd ..
Rscript ihs.r $dir $species

