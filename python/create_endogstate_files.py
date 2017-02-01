from __future__ import division

import sys

import pandas as p
import numpy as np

sys.path.append('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/') 

simulations = np.loadtxt('/msu/scratch3/m1cjg01/aer_revision_ed/final_code/pmax_smoothed.txt')
nsave=1000

results = p.DataFrame()


from progressbar import ProgressBar
pb = ProgressBar()

for i in pb(np.arange(nsave)):

    resi = p.DataFrame(simulations[i::nsave].T,  
                       columns=['endogstate_{}'.format(j+1) for j in range(28)])
    resi['sim'] = i

    results = results.append(resi)

#results = results.groupby(results.index).mean()
outdir = '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/low-me/smoothed-states/'

for col in ['endogstate_{}'.format(j+1) for j in range(28)]:
    dat = results.groupby([results.sim, results.index])[col].mean().unstack()

    np.savetxt(outdir+col+'.txt', dat.values)