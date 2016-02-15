## Copied from changFig1Maker.py on 11 Feb 2016
## Run with python/2.6.5
## Run with ssh -X
## Script calculate allele frequency spectrum
## Script presumes that there is no missing data when figuring out how many projections to use


##################################
######## IMPORT LIBRARIES ########
##################################

import sys
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg')
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg/dadi')

import dadi
import re
import argparse
import csv

##################################
############ PARSE CMD ###########
##################################

parser = argparse.ArgumentParser()
parser.add_argument("--dadi_syn_file", help="path to dadi synonymous mutations input file")
parser.add_argument("--dadi_nonsyn_file", help="path to dadi nonsyn mutations input file")
parser.add_argument("--dadi_genic_file", help="path to dadi genic mutations input file")
parser.add_argument("--dadi_intergenic_file", help="path to dadi intergenic mutations input file")
parser.add_argument("--ms_file", help="path to ms simulation results file")
parser.add_argument("--out", help="define an output AFS file name")
args = parser.parse_args()

##################################
######## PRODUCE SPECTRUM ########
##################################

## determine number of samples (so I can determine number of projections needed)
import csv
file = open(args.dadi_syn_file, 'r')
tsv = csv.reader(file, delimiter="\t")
next(tsv)
line = next(tsv)
num_proj = int(line[3]) + int(line[5])
	# calculate number of samples / projections
	# need to convert strings into intergers

## make data dictionaries
dd_syn = dadi.Misc.make_data_dict(args.dadi_syn_file)
dd_nsy = dadi.Misc.make_data_dict(args.dadi_nonsyn_file)
dd_gen = dadi.Misc.make_data_dict(args.dadi_genic_file)
dd_ige = dadi.Misc.make_data_dict(args.dadi_intergenic_file)

## make folded spectra
fs_syn = dadi.Spectrum.from_data_dict(dd_syn, pop_ids = ['Pop1'], projections = [num_proj], polarized=False)
fs_nsy = dadi.Spectrum.from_data_dict(dd_nsy, pop_ids = ['Pop1'], projections = [num_proj], polarized=False)
fs_gen = dadi.Spectrum.from_data_dict(dd_gen, pop_ids = ['Pop1'], projections = [num_proj], polarized=False)
fs_ige = dadi.Spectrum.from_data_dict(dd_ige, pop_ids = ['Pop1'], projections = [num_proj], polarized=False)
	# since none of my variant records have missing calls for some individuals
	# (since I dealt with that in filtering the VCF in GATK), no need to project down.
	# Projecting down can give you more SNPs if you have missing calls in some indivs


#####################################################
####### MS COAL SIMULATIONS MUST BE RUN BY NOW ######
#####################################################

## too complicated
#import os
#core = " -n 1 -n 2 - ej 0.3 2 1 "
#command = dadi.Misc.ms_command(theta=1000, ns=(20 ,20), core=core, iter=1000, recomb=0.3)
#ms_fs = dadi.Spectrum.from_ms_file(os.popen(command))

# read in the data
ms = dadi.Spectrum.from_ms_file(args.ms_file) # CP2 theta
#ms_78_1000_t368 = dadi.Spectrum.from_ms_file('dadi/ms/ms_78_1000_t368') # CP2 theta
#ms_78_10M_t000478 = dadi.Spectrum.from_ms_file('dadi/ms/ms_78_10M_t0.00478') # for pf
#ms_78_1000_s7000 = dadi.Spectrum.from_ms_file('dadi/ms/ms_78_1000_s7000') # ~ number of Pf snps
#ms_78_1000_s15000 = dadi.Spectrum.from_ms_file('dadi/ms/ms_78_1000_s15000') # not sure where this theta came from
#ms_70_1000_t11100 = dadi.Spectrum.from_ms_file('dadi/ms/ms_70_1000_t11100') # theta reflective of Pv
#ms_70_1000_s80000 = dadi.Spectrum.from_ms_file('dadi/ms/ms_70_1000_s80000') # ~ number of Pv snps



# fold it
ms_folded = ms.fold()
#ms_78_1000_t368_folded = ms_78_1000_t368.fold()
#ms_78_10M_t000478_folded = ms_78_10M_t000478.fold()
#ms_78_1000_s15000_folded = ms_78_1000_s15000.fold()
#ms_78_1000_s7000_folded = ms_78_1000_s7000.fold()
#ms_70_1000_t11100_folded = ms_70_1000_t11100.fold()
#ms_70_1000_s80000_folded = ms_70_1000_s80000.fold()

#dadi.Plotting.plot_1d_fs(ms_fs_folded, fig_num=None)
#pylab.show()


##################################
######## PRINTING THE AFS ########
##################################

AFS_OUT = open(args.out,'w')

## Add a header
print >>AFS_OUT, "syn", "\t", "nonsyn", "\t",	"gen", "\t", "ige", "\t", "ms"

for i in range(len(fs_syn)):
	print >>AFS_OUT, fs_syn[i]/fs_syn.S(), "\t", fs_nsy[i]/fs_nsy.S(), "\t", fs_gen[i]/fs_gen.S(), "\t", fs_ige[i]/fs_ige.S(), "\t", ms_folded[i]/ms_folded.S()

AFS_OUT.close()

