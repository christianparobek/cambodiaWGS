## I need to merge Miotto's ENA/EBI genomes
## The accession names are arbitrary, so I'll need a list as input
## Actually need to do two rounds of merging.
##	Merging runs down into a single library
##	Merging libraries down into a single sample
## Started 17 September 2015


#########################################
##### DECLARE SOME USEFUL VARIABLES #####
#########################################

# metadata directory
meta='/proj/julianog/users/ChristianP/cambodiaWGS/extraPops/Miotto2013'

###########################################################
##### FIRST MERGE MULTIPLE RUNS INTO A SINGLE LIBRARY #####
###########################################################

###########################################################
#### BECAUSE 47 DUPLICATE ASIAN LIBS ARE MISSING FROM #####
#### MIOTTO'S METADATA THIS WILL MERGE TO SAMPLE LEVEL ####
###########################################################

# declare arrays to keep the ers (libs) and corresponding errs (runs)
ers=()
err=()

for ERS in `cat $meta/miottoAsiaMetadata.txt | grep "ERS" | cut -f2 | sort | uniq`
do

	# push ers names to array
	ers+=($ERS)

	# push err names to array
	err+=(`grep $ERS $meta/miottoAsiaMetadata.txt | cut -f7 | tr '\n' ',' | sed 's/,$//'`)

done


# iterate through each ers/errs pair and either do symlink or merge
# these "runmerge" files are all runs merged into a single library
for i in $(eval echo "{1..${#ers[*]}}")
do

	# split errs from one library into their own array
	errs=(`echo ${err[$i-1]} | tr ',' '\n'`)
	echo ${#errs[*]}

	# setup variables for picard merge command
	# to be used in elif expression below
	beginning='java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MergeSamFiles.jar'
	ending='SORT_ORDER=coordinate MERGE_SEQUENCE_DICTIONARIES=true'
	ins=''
	out="O=aln/${ers[$i-1]}.merged.bam"

	# create symlink for all libraries with only one run
	if test ${#errs[*]} = 1
	then
		abspath=$(readlink -e aln/`ls aln | grep ${err[$i-1]} | grep dedup.bam`)
		ln -s $abspath "aln/${ers[$i-1]}.merged.bam"

	# merge all libs with multiple runs
	elif test ${#errs[*]} > 1
	then
		# for each err (run) in this ers (library), make an input term
		for i in $(eval echo "{1..${#errs[*]}}")
		do
			ins="$ins I=aln/${errs[$i-1]}_2013-04-28.dedup.bam"
		done
		# echo the merge expression
		beginning="$beginning $ins $out $ending"
		bsub $beginning

	fi
done
