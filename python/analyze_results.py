from __future__ import division
import numpy as np
import pandas as p
import sys
sys.path.append('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/')

from simulations import glss_names, from_matt
from model import texnames

import datetime as dt

started = dt.datetime(2016, 6, 13, hour=16, minute=34)
now = dt.datetime.now()


# old parameters
sim_dir =  '/msu/scratch3/m1cjg01/aer_revision_ed/final-results/'
parasim = p.read_csv(sim_dir+'simple_mcmc_results_agg.txt', delim_whitespace=True, 
                     names=glss_names+['posterior'], index_col=False)
parasim = parasim.ix[5000:]
parasim_old = from_matt(parasim) 



sim_dir = '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/low-me/'
parasim = p.read_csv(sim_dir+'parasim.txt', delim_whitespace=True, 
                     names=glss_names, index_col=False)
parasim = from_matt(parasim)
parasim = parasim.ix[1000:]
ncomplete = parasim.shape[0]
ntotal = 50000

frac = ncomplete/ntotal
elapsed = now-started
print "Simulation has ran for     :", elapsed#%i days, %i hours, and %minutes",  (elapsed.days, elapsed.hour, elapsed.minute)
print "Completed %i of 50000 draws." %  ncomplete
print "Expected Completion time   :", dt.timedelta(seconds=elapsed.total_seconds()/frac) + started
print "Average Acceptance Rate    :", parasim.beta_tr.unique().size/ncomplete

pnames = ['beta_tr', 'pibar_tr', 'gz_tr', 'alpha', 
          'gam_rs', 'gam_dp', 'gamdy', 'gamxhp', 
          'gamma', 'sigmal', 'sigmaa', 'phii', 
          'phi', 'ap_tr', 'phiw', 'aw_tr', 'rhog', 
          'sdevg_tr', 'rhoinv', 'sdevinv_tr', 
          'sdevliq_tr', 'sdevtech_tr', 'sdevint_tr', 
          'rhoelast', 'sdevelast_tr']

header = r"""
  \begin{table}[h]
    \begin{center}
      \caption{Posterior Distribution -- Linear Model}
      \vspace{0.1in}
      \begin{tabular}{lcccccc}
        \hline\hline
         & Current Baseline & Low ME + Markup  \\
        Parameter & Mean & [05, 95] & Mean & [05, 95] \\"""

footer = r"""\hline
      \end{tabular}
    \end{center}
  \end{table}"""
row =  '{:20s} & {:6.2f} & [{:6.2f}, {:6.2f}] & {:6.2f} & [{:6.2f}, {:6.2f}]  \\\\'.format

print header
for pa in pnames:
    print row(texnames[pa], parasim_old[pa].mean(), parasim_old[pa].quantile(0.05), parasim_old[pa].quantile(0.95), parasim[pa].mean(), parasim[pa].quantile(0.05), parasim[pa].quantile(0.95))
print footer


