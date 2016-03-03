## A control script to calculate and graph EHH
## for selected
## Started 2016-03-01
## USAGE: bash ehhRunner.sh <locusID> <VCF> <MAP> <outfile>


##################################
####### PARSE COMMAND LINE #######
##################################

locus=$1
vcf=$2
map=$3
out=$4


##################################
######### CALCULATE EHH ##########
##################################

selscan --ehh $locus --vcf $vcf --map $map --out $out


##################################
############ PLOT IT #############
##################################

Rscript selscan/ehh/ehh.r $out $locus 
