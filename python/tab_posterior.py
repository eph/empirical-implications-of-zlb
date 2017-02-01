import pandas as p
from itertools import izip_longest

import sys
sys.path.append('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/')
from model import texnames, mod

from simulations import from_matt, glss_names

sim_dir =  '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/final-final/'
#sim_dir = '/mq/home/m1eph00/tmp/aer_revision_ed/final_code/january-replication/'


paraest = p.DataFrame()
for i in range(1, 4):
    parasim = p.read_csv(sim_dir + 'parasim_final_{:02d}.txt'.format(i), delim_whitespace=True, 
                         names=glss_names, index_col=False)
    parasim = parasim.ix[int(0.2*parasim.shape[0]):, :]
    parasim['trial'] = i
    paraest = paraest.append(parasim)

#f = '/mq/home/m1eph00/tmp/aer_revision_ed/low-me-results/simple_mcmc_final_estimates.txt'
# parasim = p.read_csv(f, delim_whitespace=True, 
#                      names=glss_names+['post'], index_col=False)
# parasim = parasim.ix[int(0.2*parasim.shape[0]):, :]
# parasim['trial'] = 0
# paraest = paraest.append(parasim)

paraest = from_matt(paraest) 

header = r"""
  \begin{table}[h]
    \begin{center}
      \caption{Posterior Distribution -- Nonlinear Model}
      \vspace{0.1in}
      \begin{tabular}{lcclcc}
        \hline\hline
        Parameter & Mean & [05, 95] &         Parameter & Mean & [05, 95] \\
"""

footer = r"""\hline
      \end{tabular}
    \end{center}
  \end{table}"""

row = r"""
        \hline
        \multicolumn{{6}}{{c}}{{{}}} \\
        \hline""".format

def tab_section(paralist):
    out = '{:20s} & {:5.2f} ({:5.2f}) & [{:5.2f},  {:5.2f}] '.format
    text = lambda x: out(texnames[x], paraest[paraest.trial < 3][x].mean(), 
                         paraest[x].groupby(paraest.trial).mean().std(), 
                         *paraest[paraest.trial<3][x].quantile([0.05, 0.95]))
    # out = '{:20s} & {:5.2f} & [{:5.2f},  {:5.2f}] '.format
    # text = lambda x: out(texnames[x], paraest[x].mean(), *paraest[x].quantile([0.05, 0.95]))
    #
    paralist = map(text, paralist)
    rows = [' & '.join(r) + '\\\\' for r in izip_longest(*[paralist[i::2] for i in range(2)], 
                                                fillvalue = '& &')]
    return '\n'.join(rows)


print header
print row('Steady State')
print tab_section(['beta_tr', 'pibar_tr', 'gz_tr', 'alpha'])
print row('Policy Rule')
print tab_section(['gam_rs', 'gam_dp', 'gamdy', 'gamxhp'])
print row('Endogenous Propagation')
print tab_section(['gamma', 'sigmal', 'sigmaa', 'phii', 
                   'phi', 'ap_tr', 'phiw', 'aw_tr'])
print row('Exogenous Processes')
print tab_section(['rhog', 'sdevg_tr', 'rhoinv', 'sdevinv_tr', 
                   'sdevliq_tr', 'sdevtech_tr', 'sdevint_tr'])
print footer




