## Likely have to realign for LUMPY variant detection
## SVTyper is having major trouble with the headers of my
## PICARD & GATK pipeline BAMs
## So realign all using bwa mem and samtools

workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/pf/'
REF = '/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'
readWD = '/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/pf/'
SAMPLES, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/pf/catreads/{sample}_R1.fq.gz')



####### Target #######
rule all:
	#input: expand('aln/{sample}.bam', sample = SAMPLES)
	input: expand('aln/{sample}.discordants.bam', sample = SAMPLES)
	#input: expand('aln/{sample}.sorted.bam.bai', sample = SAMPLES)
	#input: expand('aln/{sample}.splitters.bam', sample = SAMPLES)
	#input: expand('sv/{sample}.gt.vcf', sample = SAMPLES)

rule svtyper:
	input: bam = 'aln/{sample}.sorted.bam', splitters = 'aln/{sample}.splitters.sorted.bam', vcf = 'sv/{sample}.vcf'
	output: 'sv/{sample}.gt.vcf'
	shell: 'svtyper -B {input.bam} -S {input.splitters} -i {input.vcf} -M > {output}'

rule extract_splitreads:
	input: 'aln/{sample}.sorted.bam'
	output: 'aln/{sample}.splitters.bam'
	shell: 'module remove python; \
		addpython2.7.6; \
		samtools view -h {input} \
			| /proj/julianog/src/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
			| samtools view -Sb - \
			> {output}'

rule extract_discordants:
	input: 'aln/{sample}.sorted.bam'
	output: 'aln/{sample}.discordants.bam'
	shell: 'samtools view -b -F 1294 {input} > {output}'

rule index_bam:
	input: 'aln/{sample}.sorted.bam'
	output: 'aln/{sample}.sorted.bam.bai'
	shell: 'samtools index {input}'

rule sort_bam:
	input: 'aln/{sample}.bam'
	output: 'aln/{sample}.sorted.bam'
	shell: 'samtools sort {input} aln/{wildcards.sample}.sorted'

rule fastq_to_bam:
	input: 'catreads/{sample}_R1.fq.gz', 'catreads/{sample}_R2.fq.gz'
	output: 'aln/{sample}.bam'
	shell: 'bwa mem {REF} {readWD}{input[0]} {readWD}{input[1]} \
		-M -t 1 -v 2 -A 2 -L 15 -U 9 -T 75 \
		-k 19 -w 100 -d 100 -r 1.5 -c 10000 \
		-B 4 -O 6 -E 1 -R "@RG\tID:id\tSM:sample\tLB:lib"|\
		samtools view -Sb - > {output}'
