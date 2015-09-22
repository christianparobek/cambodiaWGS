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

## import pv vcf
vcf_reader = vcf.Reader(open('/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/our_goods_UG.pass.vcf.gz', 'r'))

## import pv chrs
pvDict = SeqIO.to_dict(SeqIO.parse("/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta", "fasta"))

## out file
out_handle = open("your_out_file.dadi", "w")

##################################
########### PARSE VCF ############
##################################

# Declare some lists that will be used later on
sites = []
chrNames = []
wt = []
mut = []
a1_count = []
a2_count = []

#new_vcf_reader = vcf_reader.fetch('Pv_Sal1_chr01', 0, 10000)
for record in vcf_reader:
#	print record
	
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
####### EXTRACT SEQUENCE #########
##################################

## Test only
#sites = [3,4,5,6,7,8,9]
#chrNames = ["Pv_Sal1_chr02","Pv_Sal1_chr02","Pv_Sal1_chr02","Pv_Sal1_chr02","Pv_Sal1_chr02","Pv_Sal1_chr02"]


out_handle = open("your_out_file.dadi", "w")

# Define Variables
sp = 'Pvivax'
og = 'Outgroup'
a1 = 'Allele1'
p1 = "Pop1"
a2 = 'Allele2'
ch = "Chrom"
po = "Position"
og_value = ''-'' ## NEED THIS TO OUTPUT '-' IN THE FINAL FILE!! WILL THIS WORK????

# print the header
out_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sp, og, a1, p1, a2, p1, ch, po))

## Match chr name in dictionary
## Extract base and -1/+1 base
for idx, val in enumerate(chrNames):
	out_handle.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\n" % (pvDict[val].seq[sites[idx]-2:sites[idx]+1], og_value, wt[idx], a1_count[idx], mut[idx], a2_count[idx], chrNames[idx], sites[idx]))


out_handle.close()

