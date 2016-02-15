#!/bin/bash

## Started 15 February 2016
## This is a launcher script for afsPrinter.py
## To take dadi input files and calculate an afs
## Christian Parobek
## Usage: bash afsPrinter.sh </path/to/vcf> <a_new_wk_dir_name> <species_name_either pv or pf>


############################################################
#################### PARSE CMD LINE ########################
############################################################

annovcf=$1
syn=$2
nonsyn=$3
genic=$4
intergenic=$5
out=$6


############################################################
############ RUN APPROPRIATE MS SIMULATIONS ################
############################################################

# Get the number of individuals in this VCF
num_indiv=`grep "CHROM" $annovcf | tr '\t' '\n'| tail -n +10 | wc -l`
# Get the rough number of SNPs in this subset of samples
num_snps=`grep "chr\|_v3" $annovcf | wc -l`

echo $num_indiv
echo $num_snps

# Run the ms coalescent simulator
ms $num_indiv 1000 -s $num_snps > dadi/ms/ms_$num_indiv\_1000_s$num_snps

############################################################
############## SET UP PYTHON ENVIRONMENT ###################
############################################################

## turns out "module" is a system function and it doesn't play well in a script
## took a while to get this figured out, but this seems to work
## need to empty whatever python module is loaded originally and replace with v2.6.5

eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash remove python`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash load python`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list` # sanity check for python version


############################################################
#################### CALCULATE AFS #########################
############################################################

python dadi/afsPrinter.py --dadi_syn_file $syn --dadi_nonsyn_file $nonsyn --dadi_genic_file $genic --dadi_intergenic_file $intergenic --ms_file dadi/ms/ms_$num_indiv\_1000_s$num_snps --out $out

