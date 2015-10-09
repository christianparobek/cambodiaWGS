## Use this scripts to convert to dadi format and fit the models
## 04 October 2015
## Christian Parobek

###################################
##### DEFINE USEFUL VARIABLES #####
###################################

ref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta


###################################
###### CLEANUP LEFTOVER DATA ######
###################################

rm modelFitting/pf/2pop/*


###################################
###### RUN THE MODEL FITTING ######
###################################

for comparison in 12 # 14 24 # 13 23 34 # set up the pairwise comparisons
do

i=${comparison:0:1} # get first pop in pairwise comp
j=${comparison:1:1} # get second pop in pairwise comp

## convert vcf to dadi format
python vcf2dadi.py \
	--ref $ref \
	--vcf1 data/pf/vcf/our.pf.syn.cp$i.recode.vcf \
	--vcf2 data/pf/vcf/our.pf.syn.cp$j.recode.vcf \
	--pop1name CP$i --pop2name CP$j --out data/pf/dadi/cp$i-cp$j.dadi 

## fit the three 2-pop models to the data
for iter in {1..100}
do

bsub python 2-pop.py \
	--dadi data/pf/dadi/cp$i-cp$j.dadi \
	--pop1name CP$i --pop2name CP$j \
	--outDir modelFitting/pf/2pop/ --outName cp$i-cp$j

done

done
