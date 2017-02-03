import pandas as p
from itertools import izip_longest

from model import texnames, mod
from simulations_aer2 import from_matt, glss_names

sim_dir =  'results/mcmc/'

paraest = p.DataFrame()
for i in range(5):
    parasim = p.read_csv(sim_dir + 'output{:i}_para_tmp.txt'.format(i), delim_whitespace=True, 
                         names=glss_names, index_col=False)
    parasim = parasim.ix[int(0.2*parasim.shape[0]):, :]
    parasim['trial'] = i
    paraest = paraest.append(parasim)

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
    # if you want NSE
    # out = '{:20s} & {:5.2f} ({:5.2f}) & [{:5.2f},  {:5.2f}] '.format
    # text = lambda x: out(texnames[x], paraest[x].mean(), 
    #                      paraest[x].groupby(paraest.trial).mean().std(), 
    #                      *paraest[x].quantile([0.05, 0.95]))
    out = '{:20s} & {:5.2f} & [{:5.2f},  {:5.2f}] '.format
    text = lambda x: out(texnames[x], paraest[x].mean(), *paraest[x].quantile([0.05, 0.95]))
    
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




