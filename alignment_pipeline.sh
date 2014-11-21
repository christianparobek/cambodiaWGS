### This is a Bowtie2 alignment pipeline
### Modified for Pv whole-genome sequencing
### Started August, 2014
### Features:
###	Paired-end alignments using BWA-MEM    
###	Paired-end alignments using Bowtie2
###	Variant-calling using GATK


##########################################################################
###################### REFERENCE GENOME PREPARATION ######################
##########################################################################

## INDEX REFERENCE SEQUENCE FOR BWA
#bwa 

## INDEX REFERENCE SEQUENCE FOR BOWTIE2
#bowtie2-build PvSal1_v10.0.fasta PvSal1_10.0

## INDEX REFERENCE SEQUENCE FOR SAMTOOLS... necessary for the mpileup step
#samtools faidx PvSal1_v10.0.fasta

## INDEX REFERENCE SEQUENCE FOR GATK
#

##########################################################################
###################### SAMPLE ALIGNMENT & CLEANING #######################
##########################################################################

# Aligning so many files in an automated way will be tricky.
# Will put files from all runs (lanes) into a single dir.

ref=/proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta
picard=/nas02/apps/picard-1.88/picard-tools-1.88
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar

#for name in `cat samplenames6.txt`
#do

### ALIGN PAIRED-END READS WITH BWA_MEM
#	for rgid in `cat rgidnames.txt`
#	do
#		if test -f reads/$name$rgid\_R1.fastq.gz
#		then
#			bwa mem -M \
#			-t 8 \
#			-v 2 \
#			-R "@RG\tID:$name$rgid\tPL:illumina\tLB:$name\tSM:$name" \
#			$ref \
#			reads/$name$rgid\_R1.fastq.gz \
#			reads/$name$rgid\_R2.fastq.gz \
#			> aln/$name$rgid.sam
#				# -M marks shorter split hits as secondary
#				# -t indicates number of threads
#				# -v 2 is verbosity ... warnings and errors only
#		fi
#	done

### MERGE, SORT, AND COMPRESS SAM FILES
###	Need to merge the correct number of files for each sample
###	Construct conditionals to test number of SAM files, then merge

#array=(`ls aln/ | grep $name`)

#	## One sample run on two lanes 
#	if test ${#array[*]} = 2
#	then
#		java -jar $picard/MergeSamFiles.jar \
#		I=aln/${array[0]} \
#		I=aln/${array[1]} \
#		O=aln/$name.merged.bam \
#		SORT_ORDER=coordinate \
#		MERGE_SEQUENCE_DICTIONARIES=true
#	fi

#	## One sample run on four lanes 
#	if test ${#array[*]} = 4
#	then
#		java -jar $picard/MergeSamFiles.jar \
#		I=aln/${array[0]} \
#		I=aln/${array[1]} \
#		I=aln/${array[2]} \
#		I=aln/${array[3]} \
#		O=aln/$name.merged.bam \
#		SORT_ORDER=coordinate \
#		MERGE_SEQUENCE_DICTIONARIES=true
#	fi


## ALIGN PAIRED-END READS WITH BWA_MEM - not fancy... just simple

#for name in `cat BB-KP_names.cleanup.txt`
#do

#bwa mem -M \
#	-t 4 \
#	-v 2 \
#	-R "@RG\tID:$name\tPL:illumina\tLB:$name\tSM:$name" \
#	$ref \
#	reads/$name\_HWI-ST218_402_1_R1.fastq.gz \
#	reads/$name\_HWI-ST218_402_1_R2.fastq.gz \
#	> aln/$name.sam
#		# -M marks shorter split hits as secondary
#		# -t indicates number of threads
#		# -v 2 is verbosity ... warnings and errors only

### SORT SAM FILE AND OUTPUT AS BAM
#java -jar $picard/SortSam.jar \
#	I=aln/$name.sam \
#	O=aln/$name.sorted.bam \
#	SO=coordinate

### MARK DUPLICATES
#java -jar $picard/MarkDuplicates.jar \
#	I=aln/$name.sorted.bam \
#	O=aln/$name.dedup.bam \
#	METRICS_FILE=aln/$name.dedup.metrics \
#	REMOVE_DUPLICATES=False

### INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar $picard/BuildBamIndex.jar \
#	INPUT=aln/$name.dedup.bam

### IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar $gatk \
#	-T RealignerTargetCreator \
#	-R $ref \
#	-L gatk.intervals \
#	-I aln/$name.dedup.bam \
#	-o aln/$name.realigner.intervals \
#	-nt 4

### PERFORM THE ACTUAL REALIGNMENT
#java -jar $gatk 
#	-T IndelRealigner \
#	-R $ref \
#	-L gatk.intervals \
#	-I aln/$name.dedup.bam \
#	-targetIntervals \
#	aln/$name.realigner.intervals \
#	-o aln/$name.realn.bam

#done

##########################################################################
############################ VARIANT CALLING #############################
##########################################################################

### MULTIPLE-SAMPLE VARIANT CALLING (UG IS GATK'S CALLER FOR NON-DIPLOID)
#java -jar $gatk \
#	-T UnifiedGenotyper \
#	-R $ref \
#	-L gatk.intervals \
#	-I bamnames.list \
#	-o variants/combined.vcf \
#	-ploidy 1 \
#	-nt 8

##########################################################################
########################### VARIANT FILTERING ############################
##########################################################################

### FILTER BY DEPTH IN PERCENTAGE OF SAMPLES
#java -Xmx2g -jar $gatk \
#	-T CoveredByNSamplesSites \
#	-R $ref \
#	-V variants/combined.vcf \
#	-out 10xAT80%.intervals \
#	-minCov 10 \
#	-percentage 0.80 \
#		# Output interval file contains sites that passed

## FILTER VCF BY NEAFSEY PARALOGS, TANDEM REPEAT REGIONS, AND SPECIFIC QUALITY SCORES
java -jar $gatk \
	-T VariantFiltration \
	-R $ref \
	-V variants/combined.vcf \
	-L 10xAT80%.intervals \
	-XL neafseyExclude.intervals \
	-XL trfExclude.intervals \
	--filterExpression "QD < 5.0" \
	--filterName "QD" \
	--filterExpression "MQ < 60.0" \
	--filterName "MQ" \
	--filterExpression "FS > 10.0" \
	--filterName "FS" \
	--filterExpression "MQRankSum < -5.0" \
	--filterName "MQRankSum" \
	--filterExpression "ReadPosRankSum < -5.0" \
	--filterName "ReadPosRankSum" \
	--logging_level ERROR \
	-o variants/combined.qual.vcf
		# --logging_level ERROR suppresses any unwanted messages


##########################################################################
############################## EXTRA TOOLS ###############################
##########################################################################

## CALCULATE COVERAGE
#bedtools genomecov -ibam aln/$name.realn.bam -max 10 | grep genome > coverage/$name.cov10

## GATK DEPTH OF COVERAGE CALCUALTOR
#java -Xmx10g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
#	-T DepthOfCoverage \
#	-R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta \
#	-I bamnames.list \
#	-o coverage/allExons.cov \
#	-geneList coverage/PvSal1_v10.0_exons.refseq \
#	-L coverage/PvSal1_v10.0_exons.intervals \
#	-omitBaseOutput \
#	--minMappingQuality 20 \
#	--minBaseQuality 20
#		#-omitBaseOutput gets rid of the large by-sample-by-base output file
#		# Apparently, we can provide a refseq file of features in the genome for site-by-site analysis
#		# http://gatkforums.broadinstitute.org/discussion/1329/using-refseq-data

## COUNT READS
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T CountReads -R $ref -I aln/$name.merged.bam -rf MappingQualityZero

## SORT SAM FILE AND OUTPUT AS BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=aln/$name-lane1.sam O=aln/$name-lane1.sorted.bam SO=coordinate
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=aln/$name-lane2.sam O=aln/$name-lane2.sorted.bam SO=coordinate

## INDEX BAM FILE
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=aln/1737Pv.sorted.bam

## VALIDATE VCF FORMAT FOR GATK
#java -Xmx2g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $ref -T ValidateVariants --validationTypeToExclude ALL --variant plasmoDB_vivax_snps.vcf

## COMPARE VCF FILES
#vcftools --vcf bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bt2.vcf --diff bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bwa.vcf --out bwa_vs_bt2/compare.txt
