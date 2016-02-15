#!/bin/bash

## Started 10 February 2016
## To go from single, unannotated Pf or Pv VCF
## And make annotated VCFs, split those VCFs by
## site, then make the split VCFs into dadi format.
## This is a launcher script for vcf2dadi.py
## Christian Parobek
## Usage: bash dadiMaker.sh </path/to/vcf> <a_new_wk_dir_name> <species_name_either pv or pf>

############################################################
################# DECLARE USEFUL PATHS #####################
############################################################

# directories
#wd=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/changFig1Maker
#camWGS=/proj/julianog/users/ChristianP/cambodiaWGS

# scripts / programs
#gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar
snpEff=/nas02/home/p/r/prchrist/snpEff/snpEff.jar
SnpSift=/nas02/home/p/r/prchrist/snpEff/SnpSift.jar
vcf2dadi=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/vcf2dadi.py


############################################################
################### READ COMMAND LINE ######################
############################################################

vcf=$1
wkdir=$2
species=$3

mkdir -p dadi/data/$wkdir/vcfs
mkdir -p dadi/data/$wkdir/dadi_input
mkdir -p dadi/data/$wkdir/dadi_output_afs

############################################################
############### SET SPECIES-SPECIFIC FLAGS #################
############################################################

#### SPLIT VCF BY CHR ####

if [ $species = "pv" ]
then
	siftSpecies=PvSAL1v13
	ref=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta

elif [ $species = "pf" ]
then
	siftSpecies=Pf3D7v13
	ref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta

else
	echo "Something is wrong with your species notation. Choose either pf or pv."
fi


############################################################
################ ANNOTATE PV AND PF VCFS ###################
############################################################

# make a quick shortcut
vcfdir=dadi/data/$wkdir/vcfs/$wkdir

## Run snpEff for Pf
java -Xmx4g -jar $snpEff -v -o gatk $siftSpecies $vcf > $vcfdir.anno.vcf


############################################################
################ SPLIT VCFS BY SNP TYPE ####################
############################################################

# Get the nonsyn snps
java -jar $SnpSift filter "( EFF[*].EFFECT = 'NON_SYNONYMOUS_CODING' )" $vcfdir.anno.vcf > $vcfdir.nonsyn.vcf
# Get the syn snps
java -jar $SnpSift filter "( EFF[*].EFFECT = 'SYNONYMOUS_CODING' )" $vcfdir.anno.vcf > $vcfdir.syn.vcf
# Get the genic snps
java -jar $SnpSift filter "( EFF[*].EFFECT = 'SYNONYMOUS_CODING' ) | ( EFF[*].EFFECT = 'NON_SYNONYMOUS_CODING' )" $vcfdir.anno.vcf > $vcfdir.genic.vcf
# Get the intergenic snps
java -jar $SnpSift filter "( EFF[*].EFFECT != 'SYNONYMOUS_CODING' ) & ( EFF[*].EFFECT != 'NON_SYNONYMOUS_CODING' )" $vcfdir.anno.vcf > $vcfdir.intergenic.vcf


############################################################
################### MAKE DADI FORMAT #######################
############################################################

## turns out "module" is a system function and it doesn't play well in a script
## took a while to get this figured out, but this seems to work
## need to empty whatever python module is loaded originally and replace with v2.6.5

eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash remove python`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash load python`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list` # sanity check for python version

# make a quick shortcut
dadidir=dadi/data/$wkdir/dadi_input/$wkdir

# make dadi input files
python $vcf2dadi --ref $ref --vcf1 $vcfdir.nonsyn.vcf --out $dadidir.nonsyn.dadi
python $vcf2dadi --ref $ref --vcf1 $vcfdir.syn.vcf --out $dadidir.syn.dadi
python $vcf2dadi --ref $ref --vcf1 $vcfdir.genic.vcf --out $dadidir.genic.dadi
python $vcf2dadi --ref $ref --vcf1 $vcfdir.intergenic.vcf --out $dadidir.intergenic.dadi

