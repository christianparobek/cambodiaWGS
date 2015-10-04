## Run with python/2.6.5
## Run with ssh -X
## Script to infer demographic hx from allele frequency spectrum


##################################
######## IMPORT LIBRARIES ########
##################################

import sys
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg')
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg/dadi')
import dadi
import re


##################################
######## PRODUCE SPECTRUM ########
##################################

## make data dictionaries

dd_pf_syn = dadi.Misc.make_data_dict('pf/our.pf.syn.dadi')
dd_pf_nsy = dadi.Misc.make_data_dict('pf/our.pf.nonsyn.dadi')
dd_pf_gen = dadi.Misc.make_data_dict('pf/our.pf.genic.dadi')
dd_pf_ige = dadi.Misc.make_data_dict('pf/our.pf.intergenic.dadi')


dd_pv_syn = dadi.Misc.make_data_dict('pv/our.pv.syn.dadi')
dd_pv_nsy = dadi.Misc.make_data_dict('pv/our.pv.nonsyn.dadi')
dd_pv_gen = dadi.Misc.make_data_dict('pv/our.pv.genic.dadi')
dd_pv_ige = dadi.Misc.make_data_dict('pv/our.pv.intergenic.dadi')

## make folded spectra
fs_pf_syn = dadi.Spectrum.from_data_dict(dd_pf_syn, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pf_nsy = dadi.Spectrum.from_data_dict(dd_pf_nsy, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pf_gen = dadi.Spectrum.from_data_dict(dd_pf_gen, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pf_ige = dadi.Spectrum.from_data_dict(dd_pf_ige, pop_ids = ['Pop1'], projections = [70], polarized=False)
	# since none of my variant records have missing calls for some individuals
	# (since I dealt with that in filtering the VCF in GATK), no need to project down.
	# Projecting down can give you more SNPs if you have missing calls in some indivs

fs_pv_syn = dadi.Spectrum.from_data_dict(dd_pv_syn, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pv_nsy = dadi.Spectrum.from_data_dict(dd_pv_nsy, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pv_gen = dadi.Spectrum.from_data_dict(dd_pv_gen, pop_ids = ['Pop1'], projections = [70], polarized=False)
fs_pv_ige = dadi.Spectrum.from_data_dict(dd_pv_ige, pop_ids = ['Pop1'], projections = [70], polarized=False)



##################################
####### MS COAL SIMULATIONS ######
##################################

## too complicated
#import os
#core = " -n 1 -n 2 - ej 0.3 2 1 "
#command = dadi.Misc.ms_command(theta=1000, ns=(20 ,20), core=core, iter=1000, recomb=0.3)
#ms_fs = dadi.Spectrum.from_ms_file(os.popen(command))


## Run this command from commandline: ms 70 10000000 -t 0.00478 > ms_scratch_output
ms_70_1000_s80000 = dadi.Spectrum.from_ms_file('ms_70_1000_s80000') # for pv
ms_70_1000_s15000 = dadi.Spectrum.from_ms_file('ms_70_1000_s15000') # for pf
ms_70_10M_t000478 = dadi.Spectrum.from_ms_file('ms_70_10M_t0.00478')

ms_70_1000_s80000_folded = ms_70_1000_s80000.fold()
ms_70_1000_s15000_folded = ms_70_1000_s15000.fold()
ms_70_10M_t000478_folded = ms_70_10M_t000478.fold()

#dadi.Plotting.plot_1d_fs(ms_fs_folded, fig_num=None)
#pylab.show()


##################################
######## PRINTING THE AFS ########
##################################

PF_AFS_OUT = open('pf/pf_afs.txt','w')

## Add a header
print >>PF_AFS_OUT, "pf_syn", "\t", "pf_nonsyn", "\t",	"pf_gen", "\t", "pf_ige", "\t", "ms_70_1000_s15000_folded", "\t", "ms_70_10M_t000478_folded"	

for i in range(len(fs_pf_syn)):
	print >>PF_AFS_OUT, fs_pf_syn[i]/fs_pf_syn.S(), "\t", fs_pf_nsy[i]/fs_pf_nsy.S(), "\t", fs_pf_gen[i]/fs_pf_gen.S(), "\t", fs_pf_ige[i]/fs_pf_ige.S(), "\t", ms_70_1000_s15000_folded[i]/ms_70_1000_s15000_folded.S(), "\t", ms_70_10M_t000478_folded[i]/ms_70_10M_t000478_folded.S()

PF_AFS_OUT.close()


PV_AFS_OUT = open('pv/pv_afs.txt','w')

## Add a header
print >>PV_AFS_OUT, "pv_syn", "\t", "pv_nonsyn", "\t",	"pv_gen", "\t", "pv_ige", "\t", "ms_70_1000_s80000_folded", "\t", "ms_70_1000_s15000_folded", "\t","ms_70_10M_t000478_folded"	

for i in range(len(fs_pv_syn)):
	print >>PV_AFS_OUT, fs_pv_syn[i]/fs_pv_syn.S(), "\t", fs_pv_nsy[i]/fs_pv_nsy.S(), "\t", fs_pv_gen[i]/fs_pv_gen.S(), "\t", fs_pv_ige[i]/fs_pv_ige.S(), "\t", ms_70_1000_s80000_folded[i]/ms_70_1000_s80000_folded.S(), "\t", ms_70_1000_s15000_folded[i]/ms_70_1000_s15000_folded.S(), "\t", ms_70_10M_t000478_folded[i]/ms_70_10M_t000478_folded.S()

PV_AFS_OUT.close()
