"""
	Python script to fit one-population
	0. No Growth
	1. Expansion
	2. Contraction
	3. Bottleneck
	to observed data. 
	Argparse specifies filename where the FS is saved.
	Requires mydemos.py in the working dir to define demographic models.
"""

import sys
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.7.0-py2.6-linux-x86_64.egg')
sys.path.append('/proj/julianog/users/ChristianP/cambodiaWGS/dadi') # needed to get to mydemos

import numpy
import dadi
import argparse
from time import time
import mydemos


'''
	Parse command line
'''
parser = argparse.ArgumentParser()
parser.add_argument("--dadi", help="path to dadi file")
parser.add_argument("--pop1name", help="name for pop 1", default="Pop1")
parser.add_argument("--projection", help="project down to <> indivs", default="14", type=int)
parser.add_argument("--outDir", help="directory for results file", default="pf")
parser.add_argument("--outName", help="name for results file", default="noname")
args = parser.parse_args()

'''
	Housekeeping
'''

# for testing # dd = dadi.Misc.make_data_dict("dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi")
# for testing # fs = dadi.Spectrum.from_data_dict(dd, pop_ids = ["Pop1"], projections = [35], polarized=False)

dd = dadi.Misc.make_data_dict(args.dadi)
fs = dadi.Spectrum.from_data_dict(dd, pop_ids = [args.pop1name], projections = [args.projection], polarized=False)
nsam = fs._get_sample_sizes()
sample_sizes = [nsam[0],]
pts_l = [nsam[0]+10, nsam[0]+20, nsam[0]+30]

'''
	Bottleneck-Growth
'''

BGmod = dadi.Demographics1D.bottlegrowth
BG_ex = dadi.Numerics.make_extrap_log_func(BGmod)

upperBG = [1, 100, 5]
lowerBG = [0.001, 1, 0.01]
startBG = [0.25, 10, 0.75]
p0BG = dadi.Misc.perturb_params(startBG, fold=1, upper_bound=upperBG)

poptBG = dadi.Inference.optimize_log(p0BG, fs, BG_ex, pts_l, verbose = False, lower_bound = lowerBG, upper_bound = upperBG, maxiter = 20, epsilon = 1e-6)

BGmod = BG_ex(poptBG, sample_sizes, pts_l)
thetaBG = dadi.Inference.optimal_sfs_scaling(BGmod, fs)
BGmod *= thetaBG
llBG = dadi.Inference.ll(BGmod, fs)

## open out file
ISO_OUT = open(args.outDir + args.outName + '.bg.fit', "a")
	# for testing # ISO_OUT = open("dadi/something_special/" + "test" + '.bg.fit', "a")
ISO_OUT.write("llBG\tthetaBG\tetaD\tetaG\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\t%s\n" % (llBG, thetaBG, poptBG[0], poptBG[1], poptBG[2]))
ISO_OUT.close()

'''
   Positive Exponential Growth

nu: Ratio of contemporary to ancient population size
T: Time in the past at which growth began (in units of 2*Na
   generations)
n1: Number of samples in resulting SFS
pts: Number of grid points to use in integration.

'''

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

poptPOSG = dadi.Inference.optimize_log(p0POSG, fs, POSG_ex, pts_l, verbose = False, lower_bound = lowerPOSG, upper_bound = upperPOSG, maxiter = 20, epsilon = 1e-6)
# "popt" = optimized parameters
# A general purpose parameter optimizer.
# Seems like this is exploring the parameter space

POSGmod = POSG_ex(poptPOSG, sample_sizes, pts_l)
# This makes a simulated frequency spectrum: parameters, sample size, and array of grids
thetaPOSG = dadi.Inference.optimal_sfs_scaling(POSGmod, fs)
POSGmod *= thetaPOSG
llPOSG = dadi.Inference.ll(POSGmod, fs)
# Poisson log likelihood; first term is model and second term is data

## open out file
ISO_OUT = open(args.outDir + args.outName + '.posg.fit', "a")
	# for testing # ISO_OUT = open("dadi/something_special/" + "test" + '.posg.fit', "a")
ISO_OUT.write("llPOSG\tthetaPOSG\teta\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\n" % (llPOSG, thetaPOSG, poptPOSG[0], poptPOSG[1]))
ISO_OUT.close()

'''
	Negative Exponential Growth
'''

NEGGmod = dadi.Demographics1D.growth
NEGG_ex = dadi.Numerics.make_extrap_log_func(NEGGmod)

upperNEGG = [1,5]
lowerNEGG = [0.001, 0.01]
startNEGG = [0.5, 0.75]
p0NEGG = dadi.Misc.perturb_params(startNEGG, fold = 3, upper_bound = upperNEGG)

poptNEGG = dadi.Inference.optimize_log(p0NEGG, fs, NEGG_ex, pts_l, verbose = False, lower_bound = lowerNEGG, upper_bound = upperNEGG, maxiter = 20, epsilon = 1e-6)

NEGGmod = NEGG_ex(poptNEGG, sample_sizes, pts_l)
thetaNEGG = dadi.Inference.optimal_sfs_scaling(NEGGmod, fs)
NEGGmod *= thetaNEGG
llNEGG = dadi.Inference.ll(NEGGmod, fs)

## open out file
ISO_OUT = open(args.outDir + args.outName + '.negg.fit', "a")
	# for testing # ISO_OUT = open("dadi/something_special/" + "test" + '.negg.fit', "a")
ISO_OUT.write("llNEGG\tthetaNEGG\teta\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\n" % (llNEGG, thetaNEGG, poptNEGG[0], poptNEGG[1]))
ISO_OUT.close()

'''
	No population growth aka standard neutral model
'''

SNMmod = dadi.Demographics1D.snm
SNM_ex = dadi.Numerics.make_extrap_log_func(SNMmod)

upperSNM = [100]
lowerSNM = [1]
startSNM = [10]
p0SNM = dadi.Misc.perturb_params(startSNM, fold=1, upper_bound=upperSNM)

poptSNM = dadi.Inference.optimize_log(p0SNM, fs, SNM_ex, pts_l, verbose = False, lower_bound = lowerSNM, upper_bound = upperSNM, maxiter = 20, epsilon = 1e-6)

SNMmod = SNM_ex(poptSNM, sample_sizes, pts_l)
thetaSNM = dadi.Inference.optimal_sfs_scaling(SNMmod, fs)
SNMmod *= thetaSNM
llSNM = dadi.Inference.ll(SNMmod, fs)

## open out file
ISO_OUT = open(args.outDir + args.outName + '.snm.fit', "a")
	# for testing # ISO_OUT = open("dadi/something_special/" + "test" + '.snm.fit', "a")
ISO_OUT.write("llNEGG\tthetaNEGG\n")
ISO_OUT.write("%s\t%s\n" % (llSNM, thetaSNM))
ISO_OUT.close()

''' Positive Two Epoch Growth--Instantaneous size change some time ago
nu: Ratio of contemporary to ancient population size
T: Time in the past at which size change happened (in units of 2*Na
   generations)
n1: Number of samples in resulting SFS
pts: Number of grid points to use in integration
NFB edited 8 March 2016
'''

####Essentially copied and pasted from positive growth with the renaming of the var for Two Epoch Model
####My sense is that these initial parameters should all stay the same. Difference is in growth function/formula not population size/SFS?

EPOCHGmod = dadi.Demographics1D.two_epoch
EPOCHG_ex = dadi.Numerics.make_extrap_log_func(EPOCHGmod)
# This extrapolates to an infinitely fine grid

upperEPOCHG = [100, 5]
# Specify the upper bound for the two parameters in the EPOCHG model, would this be any different from POSG????
# These might be etaG (growth effect size) and T (time ago)
lowerEPOCHG = [1, 0.01]
startEPOCHG = [10, 0.75]
p0EPOCHG = dadi.Misc.perturb_params(startEPOCHG, fold=1, upper_bound=upperEPOCHG)
# p0 is the initial parameter set; looks like an array

poptEPOCHG = dadi.Inference.optimize_log(p0EPOCHG, fs, EPOCHG_ex, pts_l, verbose = False, lower_bound = lowerEPOCHG, upper_bound = upperEPOCHG, maxiter = 20, epsilon = 1e-6)
# "popt" = optimized parameters
# A general purpose parameter optimizer.
# Seems like this is exploring the parameter space

EPOCHGmod = EPOCHG_ex(poptEPOCHG, sample_sizes, pts_l)
# This makes a simulated frequency spectrum: parameters, sample size, and array of grids
thetaEPOCHG = dadi.Inference.optimal_sfs_scaling(EPOCHGmod, fs)
EPOCHGmod *= thetaPOSG
llEPOCHG = dadi.Inference.ll(EPOCHGmod, fs)
# Poisson log likelihood; first term is model and second term is data

## open out file
ISO_OUT = open(args.outDir + args.outName + '.two_epoch_epochg.fit', "a")
   # for testing # ISO_OUT = open("dadi/something_special/" + "test" + '.posg.fit', "a")
ISO_OUT.write("llEPOCHG\tthetaEPOCHG\teta\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\n" % (llEPOCHG, thetaEPOCHG, poptEPOCHG[0], poptEPOCHG[1]))
ISO_OUT.close()

