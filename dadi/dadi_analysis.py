## Run with python/2.6.5
## Run with ssh -X
## Script to infer demographic hx from allele frequency spectrum


##################################
######## IMPORT LIBRARIES ########
##################################

#import sys
#sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg')
#sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg/dadi')

import sys
import dadi
import re
import pylab


##################################
####### PARSE COMMAND LINE #######
##################################

## Read in a dadi file
ALLELES = sys.argv[1]


##################################
######## PRODUCE SPECTRUM ########
##################################

## make data dictionary
dd = dadi.Misc.make_data_dict(ALLELES)
## for example
#dd = dadi.Misc.make_data_dict('changFig1Maker/pv/our.pv.syn.dadi')


## make folded spectra
fs = dadi.Spectrum.from_data_dict(dd, pop_ids = ['Pop1'], projections = [70], polarized=False)
	# since none of my variant records have missing calls for some individuals
	# (since I dealt with that in filtering the VCF in GATK), no need to project down.
	# Projecting down can give you more SNPs if you have missing calls in some indivs


##################################
######## PRODUCE SPECTRUM ########
##################################

## Fit a POSG model

nsam = fs._get_sample_sizes()
sample_sizes = [nsam[0],]
pts_l = [nsam[0]+10, nsam[0]+20, nsam[0]+30]
# Array of grid sizes - smallest grid size is slightly larger than largest nsam()

POSGmod = dadi.Demographics1D.growth
POSG_ex = dadi.Numerics.make_extrap_log_func(POSGmod)
# This extrapolates to an infinitely fine grid

upperPOSG = [100, 5]
# Specify the upper bound for the two parameters in the POSG model
# These might be etaG (growth effect size) and T (time ago)
lowerPOSG = [1, 0.01]
startPOSG = [10, 0.75]
p0POSG = dadi.Misc.perturb_params(startPOSG, fold=1, upper_bound=upperPOSG)
# p0 is the initial parameter set; looks like an array

poptPOSG = dadi.Inference.optimize_log(p0POSG, fs, POSG_ex, pts_l, verbose = True, lower_bound = lowerPOSG, upper_bound = upperPOSG, maxiter = 20, epsilon = 1e-6)
# "popt" = optimized parameters
# A general purpose parameter optimizer.
# Seems like this is exploring the parameter space

POSGmod = POSG_ex(poptPOSG, sample_sizes, pts_l)
# This makes a simulated frequency spectrum: parameters, sample size, and array of grids
thetaPOSG = dadi.Inference.optimal_sfs_scaling(POSGmod, fs)
# This fits a theta value to the sfs
# Theta is a very easy parameter to fit
POSGmod *= thetaPOSG
llPOSG = dadi.Inference.ll(POSGmod, fs)
# Poisson log likelihood; first term is model and second term is data

print llPOSG, thetaPOSG, poptPOSG[0], poptPOSG[1]

dadi.Plotting.plot_1d_fs(fs, fig_num=None)
pylab.show()

##################################
####### MS COAL SIMULATIONS ######
##################################

## too complicated
#import os
#core = " -n 1 -n 2 - ej 0.3 2 1 "
#command = dadi.Misc.ms_command(theta=1000, ns=(20 ,20), core=core, iter=1000, recomb=0.3)
#ms_fs = dadi.Spectrum.from_ms_file(os.popen(command))


### Run this command from commandline: ms 70 10000000 -t 0.00478 > ms_scratch_output
#ms_70_1000_s80000 = dadi.Spectrum.from_ms_file('ms_70_1000_s80000') # for pv
#ms_70_1000_s15000 = dadi.Spectrum.from_ms_file('ms_70_1000_s15000') # for pf
#ms_70_10M_t000478 = dadi.Spectrum.from_ms_file('ms_70_10M_t0.00478')

#ms_70_1000_s80000_folded = ms_70_1000_s80000.fold()
#ms_70_1000_s15000_folded = ms_70_1000_s15000.fold()
#ms_70_10M_t000478_folded = ms_70_10M_t000478.fold()

#dadi.Plotting.plot_1d_fs(ms_fs_folded, fig_num=None)
#pylab.show()


##################################
######## PRINTING THE AFS ########
##################################

PF_AFS_OUT = open('pf_afs.txt','w')

## Add a header
print >>PF_AFS_OUT, "pf_syn", "\t", "pf_nonsyn", "\t",	"pf_gen", "\t", "pf_ige", "\t", "ms_70_1000_s15000_folded", "\t", "ms_70_10M_t000478_folded"	

for i in range(len(fs_pf_syn)):
	print >>PF_AFS_OUT, fs_pf_syn[i]/fs_pf_syn.S(), "\t", fs_pf_nsy[i]/fs_pf_nsy.S(), "\t", fs_pf_gen[i]/fs_pf_gen.S(), "\t", fs_pf_ige[i]/fs_pf_ige.S(), "\t", ms_70_1000_s15000_folded[i]/ms_70_1000_s15000_folded.S(), "\t", ms_70_10M_t000478_folded[i]/ms_70_10M_t000478_folded.S()


PF_AFS_OUT.close()


PV_AFS_OUT = open('pv_afs.txt','w')

## Add a header
print >>PV_AFS_OUT, "pv_syn", "\t", "pv_nonsyn", "\t",	"pv_gen", "\t", "pv_ige", "\t", "ms_70_1000_s80000_folded", "\t", "ms_70_1000_s15000_folded", "\t","ms_70_10M_t000478_folded"	

for i in range(len(fs_pv_syn)):
	print >>PV_AFS_OUT, fs_pv_syn[i]/fs_pv_syn.S(), "\t", fs_pv_nsy[i]/fs_pv_nsy.S(), "\t", fs_pv_gen[i]/fs_pv_gen.S(), "\t", fs_pv_ige[i]/fs_pv_ige.S(), "\t", ms_70_1000_s80000_folded[i]/ms_70_1000_s80000_folded.S(), "\t", ms_70_1000_s15000_folded[i]/ms_70_1000_s15000_folded.S(), "\t", ms_70_10M_t000478_folded[i]/ms_70_10M_t000478_folded.S()


PV_AFS_OUT.close()


##################################
######### PLAYING AROUND #########
##################################

## how many segregating sites do we have?
fs_pf_syn.S()
fs_pf_nsy.S()
fs_pf_gen.S()
fs_pf_ige.S()
	# the sum of folded = this value

## make a file with the spectrum in it
folded.to_file('folded_spectrum.txt')

## first step in making a model is to specify num of grid points
xx = dadi.Numerics.default_grid(10000)

## Then have to specify ancestral population size?
PhiManip.phi_1D ##something??

## Also have to specify integration parameters
Integration.one_pop() # takes T and nu as parameters

test_fs = dadi . Spectrum ([0 ,100 ,20 ,10 ,1 ,0])
test_fs.to_file('test_spectrum.txt')

### plot it
import pylab
dadi.Plotting.plot_single_2d_sfs(folded, vmin=0.1)
dadi.Plotting.plot_1d_fs(ms_fs_folded, fig_num=None)
pylab.show()



dadi.Plotting.plot_1d_fs(folded, fig_num=None)








