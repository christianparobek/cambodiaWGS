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
#sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/PyVCF-0.6.7-py2.6-linux-x86_64.egg')
#sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/PyVCF-0.6.7-py2.6-linux-x86_64.egg/vcf')
#sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/Counter-1.0.0-py2.6.egg')
import Bio
from Bio import SeqIO
import vcf
import argparse


##################################
############ PARSE CMD ###########
##################################

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="path to reference genome")
parser.add_argument("--vcf1", help="path to first VCF")
parser.add_argument("--vcf2", help="path to second VCF")
parser.add_argument("--vcf3", help="path to third VCF")
parser.add_argument("--popName1", help="name for pop 1", default="Pop1")
parser.add_argument("--popName2", help="name for pop 2", default="Pop2")
parser.add_argument("--popName3", help="name for pop 3", default="Pop3")
parser.add_argument("--out", help="define an output file")
args = parser.parse_args()


##################################
########### IMPORT DATA ##########
##################################

## import reference
reference = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
	# '/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'

## open out file
OUT = open(args.out, "w")
	# 'your_out_file.dadi'

## import file for reading metadata
forMetadata = vcf.Reader(open(args.vcf1, 'r'))


###################################
### DECLARE FUNCS AND VARIABLES ###
###################################

## declare some lists that will be used later on
sites = [] # position in chromosome
chrNames = [] # chromosome name
wt = [] # reference base
mut = [] # alternate base
vcf1_a1_ct = []
vcf1_a2_ct = []
vcf2_a1_ct = []
vcf2_a2_ct = []
vcf3_a1_ct = []
vcf3_a2_ct = []

## declare some variables that will be used during printing
sp = 'PlasmodiumSpp'
og = 'Outgroup'
a1 = 'Allele1'
a2 = 'Allele2'
p1 = args.popName1
p2 = args.popName2
p3 = args.popName3
ch = 'Chrom'
po = 'Position'
og_value = '\'-\'' ## Need this ti output '-' in the final file

## define some functions
def getMetadata(vcf): ## For each line/record in VCF file
   for record in vcf:
      sites.append(record.POS)
      chrNames.append(record.CHROM)
      wt.append(record.REF)
      mut.append(record.ALT[0])
         # for some reason record.REF returns str but record.ALT returns list
         # so grab the first element of it

def getDataVCF1(vcf): ## Sum the ref/alt genotypes at each position
   for record in vcf:
      geno_count = [] # declare a counter for num ref/alt at current position
      for sample in record.samples:
         geno_count.append(int(sample['GT'])) # convert to integer
      vcf1_a1_ct.append(len(geno_count)-sum(geno_count)) # sum the refs
      vcf1_a2_ct.append(sum(geno_count)) # sum the alternates

def getDataVCF2(vcf): ## Sum the ref/alt genotypes at each position
   for record in vcf:
      geno_count = [] # declare a counter for num ref/alt at current position
      for sample in record.samples:
         geno_count.append(int(sample['GT'])) # convert to integer
      vcf2_a1_ct.append(len(geno_count)-sum(geno_count)) # sum the refs
      vcf2_a2_ct.append(sum(geno_count)) # sum the alternates

def getDataVCF3(vcf): ## Sum the ref/alt genotypes at each position
   for record in vcf:
      geno_count = [] # declare a counter for num ref/alt at current position
      for sample in record.samples:
         geno_count.append(int(sample['GT'])) # convert to integer
      vcf3_a1_ct.append(len(geno_count)-sum(geno_count)) # sum the refs
      vcf3_a2_ct.append(sum(geno_count)) # sum the alternates


##################################
##### IMPORT, PROCESS, PRINT #####
##################################

## get the metadata
getMetadata(forMetadata)

## import VCFs
if args.vcf3 and args.vcf2 and args.vcf1:
   print "importing three VCFs"
   pop1_vcf = vcf.Reader(open(args.vcf1, 'r'))
   pop2_vcf = vcf.Reader(open(args.vcf2, 'r'))
   pop3_vcf = vcf.Reader(open(args.vcf3, 'r'))
   getDataVCF1(pop1_vcf)
   getDataVCF2(pop2_vcf)
   getDataVCF3(pop3_vcf)
   ## print the header
   OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sp, og, a1, p1, p2, p3, a2, p1, p2, p3, ch, po))
   ## match chr name in dictionary
   ## extract base and -1/+1 base
   for idx, val in enumerate(chrNames):
      OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%d\n" % (reference[val].seq[sites[idx]-2:sites[idx]+1], og_value, wt[idx], vcf1_a1_ct[idx], vcf2_a1_ct[idx], vcf3_a1_ct[idx], mut[idx], vcf1_a2_ct[idx], vcf2_a2_ct[idx], vcf3_a2_ct[idx], chrNames[idx], sites[idx]))
elif args.vcf2 and args.vcf1:
   print "importing two VCFs"
   pop1_vcf = vcf.Reader(open(args.vcf1, 'r'))
   pop2_vcf = vcf.Reader(open(args.vcf2, 'r'))
   getDataVCF1(pop1_vcf)
   getDataVCF2(pop2_vcf)
   ## print the header
   OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sp, og, a1, p1, p2, a2, p1, p2, ch, po))
   ## match chr name in dictionary
   ## extract base and -1/+1 base
   for idx, val in enumerate(chrNames):
      OUT.write("%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\n" % (reference[val].seq[sites[idx]-2:sites[idx]+1], og_value, wt[idx], vcf1_a1_ct[idx], vcf2_a1_ct[idx], mut[idx], vcf1_a2_ct[idx], vcf2_a2_ct[idx], chrNames[idx], sites[idx]))
elif args.vcf1:
   print "importing one VCF"
   pop1_vcf = vcf.Reader(open(args.vcf1, 'r'))
   getDataVCF1(pop1_vcf)
   ## print the header
   OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sp, og, a1, p1, a2, p1, ch, po))
   ## match chr name in dictionary
   ## extract base and -1/+1 base
   for idx, val in enumerate(chrNames):
      OUT.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\n" % (reference[val].seq[sites[idx]-2:sites[idx]+1], og_value, wt[idx], vcf1_a1_ct[idx], mut[idx], vcf1_a2_ct[idx], chrNames[idx], sites[idx]))
else:
   print "input VCFs incorrectly specified"

OUT.close()

###################################
######## PRINT DADI FORMAT ########
###################################

