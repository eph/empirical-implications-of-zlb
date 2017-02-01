from __future__ import division
import numpy as np
import pandas as p
import sys
sys.path.append('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/')

from simulations import glss_names, from_matt
from model import texnames, mod

sys.path.append('/mq/home/m1eph00/projects/dsge-book/code/helper')
from helper import SMCResults

npart, nblocks, lam = 12000, 6, 3.0

from dsge.DSGE import DSGE

funcs = {'mean':np.mean, 
         'q05': lambda x: np.percentile(x, 5), 
         'q95': lambda x: np.percentile(x, 95)}

modnew = DSGE.read('/msu/scratch3/m1cjg01/aer_revision_ed/final_code/input/model_v5_2')
model_file = '/msu/scratch3/m1cjg01/aer_revision_ed/python/linearized_model_5obs_gamxhp_JPT_priors.yaml'
modnew = DSGE.read(model_file)

# smc = SMCResults('model_v5_2_linear-mix', npart=npart, nblocks=nblocks, lam=lam, 
#                  paranames=map(str, modnew.parameters))
smc = SMCResults('GHLS_JPT_TEST-mix', npart=npart, nblocks=nblocks, lam=lam, 
                  paranames=map(str, modnew.parameters))
parasim = smc.load_draws([6])

parasim_stats = p.DataFrame({d: parasim.apply(f) for d, f in funcs.iteritems()})


smc = SMCResults('GHLS_JPT_TEST', data_dir='/mq/scratch/m1eph00/glss/', npart=4000, nblocks=3, lam=2.1, nphi=999, paranames=map(str, mod.parameters))
paraold  = smc.load_draws(range(1))
paraold['sdevelast_tr'] = np.nan
paraold_stats = p.DataFrame({d: paraold.apply(f) for d, f in funcs.iteritems()})
pnames = ['beta_tr', 'pibar_tr', 'gz_tr', 'alpha', 
          'gam_rs', 'gam_dp', 'gamdy', 'gamxhp', 
          'gamma', 'sigmal', 'sigmaa', 'phii', 
          'phi', 'ap_tr', 'phiw', 'aw_tr', 'rhog', 
          'sdevg_tr', 'rhoinv', 'sdevinv_tr', 
          'sdevliq_tr', 'sdevtech_tr', 'sdevint_tr', 'sdevelast_tr']


#result = p.concat([parasim_stats, paraold_stats], axis=1, keys=['paranew', 'paraold'])

sumstr = '{mean:6.2f} & [{q05:6.2f}, {q95:6.2f}]'.format
dfs = [parasim_stats, paraold_stats]
row = lambda x: ['{:20s}'.format(pa)]+[sumstr(**df.ix[pa].to_dict()) for df in dfs]
print '  \\\\ \n'.join([' & '.join(row(pa)) for pa in pnames])
 