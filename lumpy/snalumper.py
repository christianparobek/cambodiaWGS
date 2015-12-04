## Likely have to realign for LUMPY variant detection
## SVTyper is having major trouble with the headers of my
## PICARD & GATK pipeline BAMs
## So realign all using bwa mem and samtools

workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/'
REF = '/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'
readWD = '/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/'
SAMPLES, = glob_wildcards('/proj/julianog/users/ChristianP/cambodiaWGS/lumpy/catreads/{sample}_R1.fq.gz')



####### Target #######
rule all:
	input: expand('aln/{sample}.bam', sample = SAMPLES)


rule fastq_to_bam:
	input: 'catreads/{sample}_R1.fq.gz', 'catreads/{sample}_R2.fq.gz'
	output: 'aln/{sample}.bam'
	shell: 'bwa mem {REF} {readWD}{input[0]} {readWD}{input[1]} \
		-M -t 1 -v 2 -A 2 -L 15 -U 9 -T 75 \
		-k 19 -w 100 -d 100 -r 1.5 -c 10000 \
		-B 4 -O 6 -E 1 -R "@RG\tID:id\tSM:sample\tLB:lib"|\
		samtools view -Sb - > {output}'
