## This is an accessory script that actually runs svTyper
## I couldn't get this to work just within the lumpy.sh script
## If I execute the svtyper code bit within the lumpy.sh
## script from the pyithon/2.7.6 environment, it should fire
## off these svtyper jobs

## 09 December 2015
## Christian Parobek


#module remove python
#addpython2.7.6

samp=$1 # read in first command-line arg

svtyper -B aln/$samp.sorted.bam -S aln/$samp.splitters.sorted.bam -i sv/$samp.vcf -M > sv/$samp.gt.vcf
	# add in -M if I use the -M command in bwa mem


