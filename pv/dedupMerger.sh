for name in `ls symlinks | grep _R1 | awk '{print substr($0, 0, length($0)-23)}' | sort -u`

do

picard=/nas02/apps/picard-1.88/picard-tools-1.88
array=(`ls aln | grep $name | grep dedup.bam`)
	## One sample run on two lanes 
	if test ${#array[*]} = 2
	then
		bsub java -jar $picard/MergeSamFiles.jar \
			I=aln/${array[0]} \
			I=aln/${array[1]} \
			O=aln/$name.merged.bam \
			SORT_ORDER=coordinate \
			MERGE_SEQUENCE_DICTIONARIES=true
	fi
	## One sample run on 1 lanes 
	if test ${#array[*]} = 1
	then
		abspath=$(readlink -e aln/`ls aln | grep $name | grep dedup.bam`)
		ln -s $abspath "aln/$name.merged.bam"
	fi

done
