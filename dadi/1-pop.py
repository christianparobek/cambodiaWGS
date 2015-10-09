"""
	Python script to fit one-population
	1. Expansion
	2. Contraction
	3. Bottleneck
	to observed data. 
	Argparse specifies filename where the FS is saved.
	Requires mydemos.py in the working dir to define demographic models.
"""

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
ISO_OUT.write("llBG\tthetaBG\tetaD\tetaG\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\t%s\n" % (llBG, thetaBG, poptBG[0], poptBG[1], poptBG[2]))
ISO_OUT.close()

'''
	Positive Growth
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
ISO_OUT.write("llPOSG\tthetaPOSG\teta\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\n" % (llPOSG, thetaPOSG, poptPOSG[0], poptPOSG[1]))
ISO_OUT.close()

'''
	Negative Growth
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
ISO_OUT.write("llNEGG\tthetaNEGG\teta\tT\n")
ISO_OUT.write("%s\t%s\t%s\t%s\n" % (llNEGG, thetaNEGG, poptNEGG[0], poptNEGG[1]))
ISO_OUT.close()

'''
	No population growth
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
ISO_OUT.write("llNEGG\tthetaNEGG\n")
ISO_OUT.write("%s\t%s\n" % (llSNM, thetaSNM))
ISO_OUT.close()

