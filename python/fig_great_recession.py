### Using in aer2 revision

import numpy as np

from simulations_aer2 import load_simulation
from figures import saved_figure, despine, figure_defaults

import pandas as p

from pandas import ExcelWriter

from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib

baseline = load_simulation().groupby(level=0).mean()['2008':]
only_liq = load_simulation(sim_dir='results/alt-sims/great_recession_liq').groupby(level=0).mean()['2008':]
only_inv = load_simulation(sim_dir='results/alt-sims/great_recession_inv').groupby(level=0).mean()['2008':]
only_tech = load_simulation(sim_dir='results/alt-sims/great_recession_tech').groupby(level=0).mean()['2008':]


to_plot = ['Output Growth', 
           'Investment Growth', 
           'Consumption Growth', 
           'Inflation', 
           'Interest Rate', 
           'Notional Rate']

col = figure_defaults()

with saved_figure('drivers_of_gr.pdf', nrows=3, ncols=2) as (fig, ax):

    for axi, obs in zip(ax.reshape(-1), to_plot):

        baseline[obs].plot(ax=axi, linewidth=2)
        only_liq[obs].plot(ax=axi, linewidth=2, linestyle='dashed')
        only_inv[obs].plot(ax=axi, linewidth=2, linestyle='dotted')
        only_tech[obs].plot(ax=axi, linewidth=2, marker='*', markevery=4)

        axi.set_title(obs,fontsize=14)         
    
        axi.get_xaxis().set_minor_locator(MultipleLocator(4))
        axi.tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
       # axi.set_xticklabels(p.PeriodIndex(range(2008,2015,1),freq='A'))
        axi.set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
        axi.set_xticklabels(['2008','2009','2010','2011','2012','2013'])
        axi.tick_params('x', length=0,width=0,which='major',direction='in',labelsize=10,pad=5)
        axi.tick_params('y',labelsize=14)
        axi.set_xlim(['2008','2014'])  
        
    ax[0,0].set_yticks(np.arange(-4,3,2))
    ax[0,0].set_ylim([-4,2])
    ax[0,0].grid(False)
        
    ax[0,1].set_yticks(np.arange(-9,7,3))
    ax[0,1].set_ylim([-9,6])
    ax[0,1].grid(False)
        
    ax[1,0].set_yticks(np.arange(-2,3,1))
    ax[1,0].set_ylim([-2,2])
    ax[1,0].grid(False)
        
    ax[1,1].set_yticks(np.arange(-1,5,1))
    ax[1,1].set_ylim([-1,4])
    ax[1,1].grid(False)
         
    ax[2,0].set_yticks(np.arange(-2,7,2))
    ax[2,0].set_ylim([-2,6])
    ax[2,0].grid(False)
        
    ax[2,1].set_yticks(np.arange(-6,7,3))
    ax[2,1].set_ylim([-6,6])
    ax[2,1].grid(False)    
    
     
    despine()     
    
    #fig.legend(,to_plot,loc=(.5,.5),ncol=1)
    l1 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0])
    l2 = matplotlib.lines.Line2D([0],[0],linestyle='dashed',linewidth=2,color=col[1])
    l3 = matplotlib.lines.Line2D([0],[0],linestyle='dotted',linewidth=2,color=col[2])
    l4 = matplotlib.lines.Line2D([0],[0],marker='*',linewidth=2,color=col[3])
    
    fig.legend([l1,l2,l3,l4],['Baseline', 'Only Risk Premium Shocks', 'Only MEI Shocks', 'Only Tech. Shocks'],bbox_to_anchor=[0.5,0.02],loc='center',ncol=2,fontsize=10) 
        

    
    
    fig.tight_layout()
    #ax[2, 0].legend(['Baseline', 'Only Liquidity Shocks', 'Only MEI Shocks', 'Only Tech. Shocks'], bbox_to_anchor=(0.5, 0.02), ncol=4)
    #fig.subplots_adjust(bottom=0.27)
    fig.subplots_adjust(hspace = .5, wspace = .4)
    fig.set_size_inches(6,8)
