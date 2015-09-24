## Script to turn vcf into dadi format
## Started 27 August 2015
## Coded in python/2.6 env
## MUST export the PYTHONPATH.... 
## export PYTHONPATH=$PYTHONPATH:/nas02/home/p/r/prchrist/Library/lib/python2.6/site-packages
## adding to sys.path.append doesn't find all the modules it needs

##################################
######## IMPORT LIBRARIES ########
##################################

## OR THIS WORKS BUT IT's CLUNKY
import sys
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/PyVCF-0.6.7-py2.6-linux-x86_64.egg')
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/PyVCF-0.6.7-py2.6-linux-x86_64.egg/vcf')
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/Counter-1.0.0-py2.6.egg')
import Bio
from Bio import SeqIO
import vcf


##################################
####### IMPORT or PARSE CMD ######
##################################

REFERENCE = sys.argv[1]
INFILE = sys.argv[2]
OUTFILE = sys.argv[3]

## import reference
pvDict = SeqIO.to_dict(SeqIO.parse(REFERENCE, "fasta"))
	# '/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'

## import VCF
vcf_reader = vcf.Reader(open(INFILE, 'r'))
	# '/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/our_goods_UG.pass.vcf.gz'

## out file
out_handle = open(OUTFILE, "w")
	# 'your_out_file.dadi'

##################################
########### PARSE VCF ############
##################################

# Declare some lists that will be used later on
sites = [] # position in chromosome
chrNames = [] # chromosome name
wt = [] # reference base
mut = [] # alternate base
a1_count = []
a2_count = []

# For each line/record in VCF file
for record in vcf_reader:
	sites.append(record.POS)
	chrNames.append(record.CHROM)
	wt.append(record.REF)
	mut.append(record.ALT[0])
		# for some reason record.REF returns str but record.ALT returns list
		# so grab the first element of it
	
	geno_count = [] # declare a counter for num ref/alt at current position
	for sample in record.samples:
		geno_count.append(int(sample['GT'])) # convert to integer
	
	a1_count.append(len(geno_count)-sum(geno_count)) # sum the refs
	a2_count.append(sum(geno_count)) # sum the alternates


##################################
####### PRINT DADI FORMAT ########
##################################

#out_handle = open("your_out_file.dadi", "w")

# Define Variables
sp = 'PlasmodiumSpp'
og = 'Outgroup'
a1 = 'Allele1'
p1 = "Pop1"
a2 = 'Allele2'
ch = "Chrom"
po = "Position"
og_value = '\'-\'' ## Need this ti output '-' in the final file

# print the header
out_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sp, og, a1, p1, a2, p1, ch, po))

## Match chr name in dictionary
## Extract base and -1/+1 base
for idx, val in enumerate(chrNames):
	out_handle.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\n" % (pvDict[val].seq[sites[idx]-2:sites[idx]+1], og_value, wt[idx], a1_count[idx], mut[idx], a2_count[idx], chrNames[idx], sites[idx]))

out_handle.close()

