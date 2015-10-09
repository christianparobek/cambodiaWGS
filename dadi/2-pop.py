"""
	Python script to fit two-population
	1. Split with exchange
	2. Split with isolation
	3. No split 
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
parser.add_argument("--pop2name", help="name for pop 2", default="Pop2")
parser.add_argument("--outDir", help="directory for results file", default="pf")
parser.add_argument("--outName", help="name for results file", default="noname")
args = parser.parse_args()

'''
	Housekeeping
'''

dd = dadi.Misc.make_data_dict(args.dadi)
fs = dadi.Spectrum.from_data_dict(dd, pop_ids = [args.pop1name, args.pop2name], projections = [14, 14], polarized=False)
nsam = fs._get_sample_sizes()
sample_sizes = [nsam[0],nsam[0]]
pts_l = [nsam[0]+10, nsam[0]+20, nsam[0]+30]


'''
	This is the SNM model - no divergence
'''
T1 = time()
SNMmod = mydemos.snm2
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
T2 = time()

## open out file
SNM_OUT = open(args.outDir + args.outName + '.snm.fit', "a")
SNM_OUT.write("llSNM\tthetaSNM\truntime\n")
SNM_OUT.write("%s\t%s\t%s\n" % (llSNM, thetaSNM, T2-T1))
SNM_OUT.close()

'''
	This is the ISO model - divergence with no gene flow
'''
T1 = time()
ISOmod = mydemos.Iso
ISO_ex = dadi.Numerics.make_extrap_log_func(ISOmod)

upperISO = [0.9999, 20, 20, 5] # originally 0.99, 20, 20, 5
lowerISO = [0.0001, 0.5, 0.5, 0.005] # originally 0.01, 0.5, 0.5, 0.005
startISO = [0.25, 5, 5, 0.75]
p0ISO = dadi.Misc.perturb_params(startISO, fold=1, upper_bound=upperISO)

poptISO = dadi.Inference.optimize_log(p0ISO, fs, ISO_ex, pts_l, verbose = False, lower_bound = lowerISO, upper_bound = upperISO, maxiter= 20, epsilon = 1e-6)

ISOmod = ISO_ex(poptISO, sample_sizes, pts_l)
thetaISO = dadi.Inference.optimal_sfs_scaling(ISOmod, fs)
ISOmod *= thetaISO
llISO = dadi.Inference.ll(ISOmod, fs)
T2 = time()

## open out file
ISO_OUT = open(args.outDir + args.outName + '.iso.fit', "a")
ISO_OUT.write("llISO\tthetaISO\ts\teta1\teta2\tT\truntime\n")
ISO_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (llISO, thetaISO, poptISO[0], poptISO[1], poptISO[2], poptISO[3], T2-T1))
ISO_OUT.close()


'''
	This is the IM model - divergence with gene flow
'''
T1 = time()
IMmod = mydemos.IM
IM_ex = dadi.Numerics.make_extrap_log_func(IMmod)

upperIM = [0.9999, 20, 20, 5, 20, 20] # originally 0.99, 20, 20, 5, 20, 20
lowerIM = [0.0001, 0.5, 0.5, 0.005, .01, .01] # originally 0.01, 0.5, 0.5, 0.005, 0.1, 0.1
startIM = [0.25, 5, 5, 0.75, 2, 2]
p0IM = dadi.Misc.perturb_params(startIM,fold=1,upper_bound=upperIM)

poptIM = dadi.Inference.optimize_log(p0IM, fs, IM_ex, pts_l, verbose = False, lower_bound = lowerIM, upper_bound = upperIM, maxiter=20, epsilon = 1e-6)

IMmod = IM_ex(poptIM, sample_sizes, pts_l)
thetaIM = dadi.Inference.optimal_sfs_scaling(IMmod, fs)
IMmod *= thetaIM
llIM = dadi.Inference.ll(IMmod, fs)
T2 = time()

## open out file
IM_OUT = open(args.outDir + args.outName + '.im.fit', "a")
IM_OUT.write("llIM\tthetaIM\ts\teta1\teta2\tT\tm12\tm21\truntime\n")
IM_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (llIM, thetaIM, poptIM[0], poptIM[1], poptIM[2], poptIM[3], poptIM[4], poptIM[5], T2-T1))
IM_OUT.close()
