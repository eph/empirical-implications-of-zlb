from model import texnames, mod
from simulations import paraest
from itertools import izip_longest
from figures import saved_figure, despine, figure_defaults


import numpy as np
import pandas as p

from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib

from pandas import ExcelWriter


sim_dir =  'results/mcmc/'

paraest = p.DataFrame()
for i in range(5):
    parasim = p.read_csv(sim_dir + 'output{:i}_para_tmp.txt'.format(i), delim_whitespace=True, 
                         names=glss_names, index_col=False)
    parasim = parasim.ix[int(0.2*parasim.shape[0]):, :]
    parasim['trial'] = i
    paraest = paraest.append(parasim)

nonlinear = from_matt(paraest) 


import sys
sys.path.append('/msu/home/m1eph00/projects/dsge-book/code/helper/')
from helper import SMCResults

smc = SMCResults('GHLS_JPT_TEST', 
                 data_dir='/msu/scratch3/m1cjg01/aer_revision_ed/python/glss/', 
                 npart=4000, nblocks=3, lam=2.1, nphi=999, paranames=map(str, mod.parameters))
linear = smc.load_draws(range(1))

col=figure_defaults()

with saved_figure('linear_vs_nonlinear.pdf', nrows=1, ncols=2) as (fig, ax):

    for axi, pa, xlim in zip(ax.reshape(-1), ['phii', 'sdevinv_tr'], [(0, 10), (0, 20)]):
        nonlinear[pa].plot(ax=axi, kind='density')
        linear[pa].plot(ax=axi, linestyle='dashed', kind='density')
        axi.set_title(texnames[pa], fontsize=16)
        axi.set_xlim(xlim)
    
    # ax[0].set_title('Output',fontsize=16)
    # #ax[0].get_xaxis().set_minor_locator(MultipleLocator(4))
    # #ax[0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
    
    ax[0].set_xticks(np.arange(0,11,2))
    ax[0].tick_params('x', length=5,width=1,which='major',direction='in',labelsize=14,pad=5)
    #ax[0].set_xticklabels([])
    ax[0].tick_params('y',labelsize=14)
    ax[0].set_xlim(0,10)
    ax[0].annotate('Linear', xy=(1,0.42),fontsize=14,color=col[1])
    ax[0].annotate('Nonlinear', xy=(5,0.2),fontsize=14,color=col[0])
    ax[0].set_ylim(0,0.5)
    ax[0].set_ylabel('')
    
    # #ax[1].get_xaxis().set_minor_locator(MultipleLocator(4))
    # #ax[1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
    ax[1].set_xticks(np.arange(0,21,5))
    ax[1].tick_params('x',length=5,width=1,which='major',direction='in',labelsize=14,pad=5)
    #ax[1].set_xticklabels([])
    ax[1].tick_params('y',labelsize=14)
    ax[1].set_xlim(0,20)
    ax[1].set_ylim(0,0.75)
    ax[1].set_ylabel('')
    
    despine()
    
    # #fig.tight_layout()
    fig.set_size_inches(8, 4)
    fig.subplots_adjust(wspace = .25)

    data_phii_nonlinear = p.DataFrame(ax[0].get_lines()[0].get_xydata())
    data_phii_linear = p.DataFrame(ax[0].get_lines()[1].get_xydata())
    data_sdevinv_nonlinear = p.DataFrame(ax[1].get_lines()[0].get_xydata())
    data_sdevinv_linear = p.DataFrame(ax[1].get_lines()[1].get_xydata())


# 508
# writer = ExcelWriter('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/linear_vs_nonlinear.xlsx')
# data_phii_nonlinear.to_excel(writer,sheet_name='phii_nonlinear')
# data_phii_linear.to_excel(writer,sheet_name='phii_linear')
# data_sdevinv_nonlinear.to_excel(writer,sheet_name='stdevinv_nonlinear')
# data_sdevinv_linear.to_excel(writer,sheet_name='stdevinv_linear')
# writer.save()
