## A Snakefile


##########################################################################################

####### Turn on for Pv ##########
workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/pv/'
REF = '/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'
readWD = '/proj/julianog/users/ChristianP/cambodiaWGS/pv/'
DATEDSAMPS, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/pv/symlinks/{ds}_R1.fastq.gz')
SAMPLES, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/pv/aln/{sample}.merged.bam')
CHROMOS, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/pv/intervals/{sample}.merged.bam')
merger = 'pv/pvDedupMerger.sh'

######## Always on #########
PICARD = '/nas02/apps/picard-1.88/picard-tools-1.88'
#GATK = '/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar'
GATK = '/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar'
TMPDIR = '/netscr/prchrist/tmp_for_picard/'


##########################################################################################


####### Target #######
rule all:
#	input: expand('aln/{ds}.dedup.bam', ds = DATEDSAMPS) # Run to here frist, then run dedupMerger.sh
#	input: 'names/bamnames.list'
#	input: 'coverage/cov_plot.pdf'
#	input: 'coverage/coverage.txt'
	input: expand('variants/{list}_UG.pass.vcf', list = 'our_goods all_goods'.split())
#	input: expand('variants/indivs/{names}.vcf', names = SAMPLES)


#rule split_vcf :
#	input: vcf = 'variants/all_goods_UG.pass.vcf', names = SAMPLES
#	output: 'variants/indivs/{names}.vcf'
#	shell: 'java -jar {GATK} \
#		-T SelectVariants \
#		-R {REF} \
#		--variant {input.vcf} \
#		-sn {input.names} \
#		-o {output}'

rule select_variants :
	input: 'variants/{list}_UG.qual.vcf'
	output: 'variants/{list}_UG.pass.vcf'
	shell: 'java -jar {GATK} \
		-T SelectVariants \
		-R {REF} -V {input} -o {output} \
		-select "vc.isNotFiltered()" \
		-restrictAllelesTo BIALLELIC'
			# keeps only unfiltered sites

rule filter_variants :
	input: vcf = 'variants/{list}_UG.vcf', depth = 'intervals/{list}_UG_05xAT100%.intervals', para = 'intervals/neafseyExclude.intervals', trf = 'intervals/trfExclude.intervals', telo = 'intervals/subtelomeres.intervals'
	output: 'variants/{list}_UG.qual.vcf'
	shell: 'java -jar {GATK} \
		-T VariantFiltration \
		-R {REF} \
		-V {input.vcf} \
		-L {input.depth} \
		-XL {input.para} \
		-XL {input.trf} \
		-XL {input.telo} \
		--filterExpression "QD < 25.0" \
		--filterName "QD" \
		--filterExpression "MQ < 55.0" \
		--filterName "MQ" \
		--filterExpression "FS > 10.0" \
		--filterName "FS" \
		--filterExpression "MQRankSum < -5.0" \
		--filterName "MQRankSum" \
		--filterExpression "ReadPosRankSum < -5.0" \
		--filterName "ReadPosRankSum" \
		--logging_level ERROR \
		-o {output}'

rule filter_min_depth :
	input: 'variants/{list}_UG.vcf'
	output: 'intervals/{list}_UG_05xAT100%.intervals'
	shell: 'java -Xmx2g -jar {GATK} \
		-T CoveredByNSamplesSites \
		-R {REF} -V {input} -out {output} \
		-minCov 05 -percentage 0.99999'
		 #Output interval file contains sites that passed
		 #Would be more elegant to use 1.0 instad of 0.99999, but that doesn't work

rule unified_genotyper :
	input: bams = 'names/{list}_5x@80%.list', intervals = 'intervals/all_chrs.intervals'
	output: 'variants/gvcfs/{list}_UG.vcf'
	shell: 'java -jar {GATK} -T UnifiedGenotyper \
		-R {REF} -I {input.bams} \
		-L {input.intervals} -nt 8 \
		-ploidy 1 -o {output}'
		# all_chrs.intervals includes only chrs and mito

rule plot_coverage:
	input: 'coverage/covPlotter.r'
	output: 'coverage/cov_plot.pdf'
	shell: 'Rscript {input}'

rule make_R_cov_script:
	input: 'coverage/coverage.txt'
	output: 'coverage/covPlotter.r'
	shell: """echo '## Read in data
		coverage <- read.table("coverage/coverage.txt", header=FALSE)

		## Order data acendingly
		coverage <- coverage[with(coverage, order(V2)), ]

		## Add a third column with numbers in it
		coverage$V4 <- 1:length(coverage$V1)

		##Plot the coverage graph
		pdf(file="coverage/cov_plot.pdf", width = 28, height = 4)
		plot(coverage$V2 ~ coverage$V4, axes=FALSE, xlab="", ylab="Frac Genome Covered", ylim=c(0,1), col="black", pch=20)
		points(coverage$V3 ~ coverage$V4, col="grey", pch=20)
		axis(1, at=1:length(coverage$V1), labels=coverage$V1, cex.axis=.7, las=3, cex = 0.5)
		axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), line=-4)
		abline(h=c(0.8, 0.75, 0.7, 0.65, 0.6), col="grey")
		legend(1,0.5, legend=c("5x Coverage", "10x Coverage"), col=c("black", "grey"), pch=20)
		dev.off()' > coverage/covPlotter.r
		"""

rule digest_coverage:
	input: expand('coverage/data/{sample}.{cov}', sample = SAMPLES, cov = 'cov05 cov10'.split())
	output: 'coverage/coverage.txt'
	shell: 'for name in `ls coverage/data/ | grep cov05 | sed "s/\.cov..//"`; \
		do \
		cov05=$(tail -1 coverage/data/$name.cov05 | cut -f 5); \
		cov10=$(tail -1 coverage/data/$name.cov10 | cut -f 5); \
		echo -e $name"\t"$cov05"\t"$cov10 >> coverage/coverage.txt; \
		done'

rule calculate_10x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov10'
	shell: 'bedtools genomecov \
		-ibam {input} -max 10 | grep genome \
		> {output}'

rule calculate_05x_cov:
	input: 'aln/{sample}.realn.bam'
	output: 'coverage/data/{sample}.cov05'
	shell: 'bedtools genomecov \
		-ibam {input} -max 5 | grep genome \
		> {output}'

rule realn_indels:
	input: bam = 'aln/{sample}.merged.bam', chrs = 'intervals/all_chrs.intervals', targets = 'aln/{sample}.realigner.intervals', 
	output: 'aln/{sample}.realn.bam'
	shell: 'java -jar {GATK} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-L {input.chrs} -targetIntervals {input.targets} \
		-o {output}' 
		# all_chrs.intervals includes just chrs and mito
		# -fixMisencodedQuals must be added for SRA data

rule find_indels:
	input: bam = 'aln/{sample}.merged.bam', index = 'aln/{sample}.merged.bai', chrs = 'intervals/all_chrs.intervals'
	output: 'aln/{sample}.realigner.intervals'
	shell: 'java -jar {GATK} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-L {input.chrs} -o {output}'
		# all_chrs.intervals includes just  chrs and mito

rule index_merged: 
	input: 'aln/{sample}.merged.bam'
	output: 'aln/{sample}.merged.bai'
	shell: 'java -jar {PICARD}/BuildBamIndex.jar INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'

rule mark_dups:
	input: 'aln/{ds}.bam'
	output:'aln/{ds}.dedup.bam','aln/{ds}.dedup.metrics'
	shell: 'java -jar {PICARD}/MarkDuplicates.jar \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=TRUE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'

rule sam_to_bam:
	input: 'aln/{ds}.sam'
	output: 'aln/{ds}.bam'
	shell: 'java -jar {PICARD}/SortSam.jar \
		I={input} O={output} \
		SO=coordinate TMP_DIR={TMPDIR}'

rule fastq_to_sam:
	input: 'symlinks/{ds}_R1.fastq.gz', 'symlinks/{ds}_R2.fastq.gz'
	output: 'aln/{ds}.sam'
	shell: 'bwa mem {REF} {readWD}{input[0]} {readWD}{input[1]} \
		-R "@RG\tID:bwa\tPL:illumina\tLB:{wildcards.ds}\tSM:{wildcards.ds[0]}{wildcards.ds[1]}{wildcards.ds[2]}{wildcards.ds[3]}{wildcards.ds[4]}{wildcards.ds[5]}{wildcards.ds[6]}{wildcards.ds[7]}{wildcards.ds[8]}{wildcards.ds[9]}" \
		-M -t 4 -v 2 -A 2 -L 15 -U 9 -T 75 \
		-k 19 -w 100 -d 100 -r 1.5 -c 10000 \
		-B 4 -O 6 -E 1 > {output}'
		# calling the @RG ID: 'bwa' because this resolves a clash with @PG ID
