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

rm modelFitting/pf/1pop/*


###################################
###### RUN THE MODEL FITTING ######
###################################

for i in 1 2 4
do

## convert vcf to dadi format
python vcf2dadi.py \
	--ref $pfref \
	--vcf1 data/pf/vcf/our.pf.syn.cp$i.recode.vcf \
	--pop1name CP$i --out data/pf/dadi/cp$i.dadi

## run the models

for rep in {1..50}
do

bsub python 1-pop.py \
	--dadi data/pf/dadi/cp$i.dadi \
	--pop1name CP$i --projection 14 \
	--outDir modelFitting/pf/1pop/ --outName cp$i

done
done
