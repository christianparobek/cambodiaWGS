## Use this scripts to convert to dadi format and fit the models
## Adapted 15 February 2015
## Christian Parobek
## Usage: bash dadi/modelFitting-1pop.sh </path/to/dadi/file> <wkdir> <outname> <species_pf_or_pv>


#################################
######## DEFINE VARIABLES #######
#################################

pfref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
pvref=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta


#################################
####### READ COMMAND LINE #######
#################################

dadi=$1
wkdir=$2
outname=$3
species=$4


#################################
######### COUNT INDIVS ##########
#################################

# pick any random line and sum up the number of individuals there
num_one=`head dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi | tail -1 | cut -f4`
num_two=`head dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi | tail -1 | cut -f6`
num_indiv=$(($num_one + $num_two)) # get the total number of individuals in this dadi file
num_proj=$(($num_indiv / 2)) # divide that number by two to get the number of projections
	# doing this because my afs is folded because I do not have ancestral data

###################################
## MKDIR & CLEANUP LEFTOVER DATA ##
###################################

mkdir -p $wkdir
rm $wkdir/*


############################################################
############## SET UP PYTHON ENVIRONMENT ###################
############################################################

## turns out "module" is a system function and it doesn't play well in a script
## took a while to get this figured out, but this seems to work
## need to empty whatever python module is loaded originally and replace with v2.6.5

eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash remove python`
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash load python/2.6.5` # has an older version of numpy for now
eval `/nas02/apps/Modules/$MODULE_VERSION/bin/modulecmd bash list` # sanity check for python version


###################################
###### RUN THE MODEL FITTING ######
###################################

## run the models
for rep in {1..100}
do

bsub python dadi/1-pop.py \
	--dadi $dadi \
	--projection $num_proj \
	--outDir $wkdir \
	--outName $outname

sleep 1 # slows down job submission to reduce chance of I/O clobbering each other

done
