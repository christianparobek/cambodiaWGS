## Use this scripts to convert to dadi format and fit the models
## 07 October 2015
## Christian Parobek

###################################
##### DEFINE USEFUL VARIABLES #####
###################################

pfref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
pvref=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta
scripts=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/robinson_sims/1-population/fitting

###################################
###### CLEANUP LEFTOVER DATA ######
###################################

rm modelFitting/pv/*


###################################
###### RUN THE MODEL FITTING ######
###################################

## convert vcf to dadi format
python vcf2dadi.py \
	--ref $pvref \
	--vcf1 data/pv/vcf/our.pv.syn.vcf \
	--pop1name Pvivax --out data/pv/dadi/our.pv.syn.dadi

## run the models

for rep in {1..100}
do

bsub python 1-pop.py \
	--dadi data/pv/dadi/our.pv.syn.dadi \
	--pop1name Pvivax --projection 70 \
	--outDir modelFitting/pv/ --outName Pvivax

done
