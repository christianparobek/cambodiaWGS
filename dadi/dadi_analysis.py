## Run with python/2.6.5
## Script to infer demographic hx from allele frequency spectrum

import sys
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg')
sys.path.append('/nas02/home/p/r/prchrist/lib/python2.6/site-packages/dadi-1.6.3-py2.6-linux-x86_64.egg/dadi')
import dadi


## make a data dictionary
dd = dadi.Misc.make_data_dict('your_out_file.dadi')

## make a folded spectrum
folded = dadi.Spectrum.from_data_dict(dd, pop_ids = ['Pop1'], projections = [70], polarized=False)

## how many segregating sites do we have
folded.S()

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


import pylab
dadi.Plotting.plot_single_2d_sfs(folded, vmin=0.1)
pylab.show()
