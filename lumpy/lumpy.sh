## A file to run lumpy structural variant analysis
## Started 01 December 2015
## Christian P

bamdir=/proj/julianog/users/ChristianP/cambodiaWGS/pv/aln
scripts=/proj/julianog/src/lumpy-sv/scripts
#samp=BB012
samp=KP017.merged
ref=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta
readdir=/proj/julianog/users/ChristianP/cambodiaWGS/pv/symlinks/

## Extract the discordant paired-end alignments
#samtools view -b -F 1294 $bamdir/$samp.merged.bam > $samp.discordants.unsorted.bam

## Extract the split-read alignments
#samtools view -h $bamdir/$samp.merged.bam \
#	| $scripts/extractSplitReads_BwaMem -i stdin \
#	| samtools view -Sb - \
#	> $samp.splitters.unsorted.bam

## Sort both alignments
#samtools sort $samp.discordants.unsorted.bam $samp.discordants
#samtools sort $samp.splitters.unsorted.bam $samp.splitters


## First, generate empirical insert size statistics on each library in the BAM file
#samtools view $bamdir/$samp.merged.bam \
#	| tail -n+100000 \
#	| $scripts/pairend_distro.py \
#	-r 101 \
#	-X 4 \
#	-N 10000 \
#	-o $samp.lib1.histo


## Run LUMPY with paired-end and split-reads
#lumpy \
#	-mw 4 \
#	-tt 0 \
#	-pe id:$samp,bam_file:$samp.discordants.bam,histo_file:$samp.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
#	-sr id:$samp,bam_file:$samp.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
#	| grep -v "AAKM" > $samp.vcf # careful with this grepping bit...

## Run SVTyper to annotate
#cat $samp.vcf \
#	| svtyper \
#		-B $bamdir/$samp.merged.bam \
#		-S $samp.splitters.bam \
#		> $samp.gt.vcf

#svtyper -B $bamdir/$samp.merged.bam -S $samp.splitters.bam -i $samp.vcf -o $samp.gt.vcf

#samtools view -H $bamdir/$samp.merged.bam | sed -e 's/PL:illumina\tPN/PN/' | sed -e 's/LB:.*VN/VN/' | sed -e 's/SM:.*-M/ -M/' | samtools reheader - $bamdir/$samp.merged.bam > $samp.reheadered.bam

##samtools view -H $bamdir/$samp.merged.bam | grep -v "@PG" | samtools reheader - $bamdir/$samp.merged.bam > $samp.reheadered.bam
#samtools index $samp.reheadered.bam


#bwa mem $ref /proj/julianog/sequence_reads/beckman_seq_backups/2014_11_10_BB_KP_WGS_Libraries/Fastq/BB012-BiooBarcode_31_CACGAT_R1.fastq.gz /proj/julianog/sequence_reads/beckman_seq_backups/2014_11_10_BB_KP_WGS_Libraries/Fastq/BB012-BiooBarcode_31_CACGAT_R2.fastq.gz \
#		-R "@RG\tID:bwa\tPL:illumina\tLB:sample" \
#		-M -t 4 -v 2 -A 2 -L 15 -U 9 -T 75 \
#		-k 19 -w 100 -d 100 -r 1.5 -c 10000 \
#		-B 4 -O 6 -E 1 > BB012.sam

#samtools sort -o BB012.bam


# Align it
bwa mem $ref $readdir/BB012_2014-11-10_R1.fastq.gz $readdir/BB012_2014-11-10_R2.fastq.gz -t 4 -R "@RG\tID:id\tSM:sample\tLB:lib"|\
	samtools view -Sb - \
	> BB012.bam

# Extract the discordant paired-end alignments
samtools view -b -F 1294 $samp.bam > $samp.discordants.unsorted.bam

# Extract the split-read alignments
module remove python
addpython2.7.6
samtools view -h $samp.bam \
	| $scripts/extractSplitReads_BwaMem -i stdin | more \
	| samtools view -Sb - \
	> $samp.splitters.unsorted.bam

# Sort both alignments
samtools sort $samp.discordants.unsorted.bam $samp.discordants
samtools sort $samp.splitters.unsorted.bam $samp.splitters


# First, generate empirical insert size statistics on each library in the BAM file
samtools view $samp.bam \
	| tail -n+100000 \
	| $scripts/pairend_distro.py \
	-r 101 \
	-X 4 \
	-N 10000 \
	-o $samp.lib1.histo


# Run LUMPY with paired-end and split-reads
module remove python
addpython3.3.3
lumpy \
	-mw 4 \
	-tt 0 \
	-pe id:$samp,bam_file:$samp.discordants.bam,histo_file:$samp.lib1.histo,mean:311,stdev:92,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
	-sr id:$samp,bam_file:$samp.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
	> $samp.vcf # careful with this grepping bit... # | grep -v "AAKM" > $samp.vcf

samtools sort $samp.bam $samp.sorted
samtools index $samp.sorted.bam

samtools index $samp.splitters.bam

module remove python
addpython2.7.6
svtyper -B $samp.sorted.bam -S $samp.splitters.bam -i $samp.vcf > $samp.gt.vcf
	# add in -M if I use the -M command in bwa mem



