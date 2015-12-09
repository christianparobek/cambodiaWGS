## A file to run lumpy structural variant analysis
## Started 01 December 2015
## Christian P

bamdir=/proj/julianog/users/ChristianP/cambodiaWGS/pv/aln
scripts=/proj/julianog/src/lumpy-sv/scripts
samp=BB012
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

#######################################
## EXTRACT THE SPLIT-READ ALIGNMENTS ##
#######################################

#module remove python
#addpython2.7.6
#for samp in `cat ../pv/names/our_goods_5x@80%.txt`
#do

#samtools view -h aln/$samp.bam \
#	| $scripts/extractSplitReads_BwaMem -i stdin \
#	| samtools view -Sb - \
#	> aln/$samp.splitters.bam

#done

#######################################
###### SORT THE ALIGNMENTS AGAIN ######
#######################################

#for samp in `cat ../pv/names/our_goods_5x@80%.txt`
#do

#echo $samp
#samtools sort aln/$samp.discordants.bam aln/$samp.discordants.sorted
#samtools sort aln/$samp.splitters.bam aln/$samp.splitters.sorted
#samtools index aln/$samp.splitters.sorted.bam

#done

#######################################
#### GET LIBRARY DISTRIBUTION INFO ####
#######################################

#for samp in `cat ../pv/names/our_goods_5x@80%.txt`
#do

#samtools view aln/$samp.bam \
#	| tail -n+100000 \
#	| $scripts/pairend_distro.py \
#	-r 101 \
#	-X 4 \
#	-N 10000 \
#	-o aln/$samp.histo > aln/$samp.libstats

#done


#######################################
############## RUN LUMPY ##############
#######################################

###module remove python
###addpython3.3.3
#for samp in `cat ../pv/names/our_goods_5x@80%.txt`
#do

#mean=`cat aln/$samp.libstats | cut -f1 | sed 's/[a-z].*://'`
#sd=`cat aln/$samp.libstats | cut -f2 | sed 's/[a-z].*://'`

#echo $samp
#echo $mean
#echo $sd

#lumpy \
#	-mw 4 \
#	-tt 0 \
#	-pe id:$samp,bam_file:aln/$samp.discordants.sorted.bam,histo_file:aln/$samp.histo,mean:$mean,stdev:$sd,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
#	-sr id:$samp,bam_file:aln/$samp.splitters.sorted.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
#	> sv/$samp.vcf # careful with this grepping bit... # | grep -v "AAKM" > $samp.vcf
#		## keep in mind that some pf were run with 125 chemistry, beginning with ones seqed in 2015
#		## dont think any Pv were 125 chemistry since they were all 2014 I think

#done

#samtools sort $samp.bam $samp.sorted
#samtools index $samp.sorted.bam


#######################################
############# RUN SVTYPER #############
#######################################



## Must execute from a python/2.7.6 environment 
#module remove python
#addpython2.7.6
for samp in `cat ../pv/names/our_goods_5x@80%.txt`
do

bsub bash svtyperStarter.sh $samp

done
