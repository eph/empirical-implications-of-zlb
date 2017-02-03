#### Used in aer2 figures
#### Generates figure with inflation, output from fig_observable_fit and technology shocks
#### LHS is baseline ME and RHS is low ME


import numpy as np

from model import data
from simulations_aer2 import load_simulation
from simulations import load_simulation_old
from figures import saved_figure, despine, figure_defaults

import pandas as p

from pandas import ExcelWriter

from matplotlib.ticker import MultipleLocator, FixedLocator
import matplotlib


sim_dir_low_me = 'results/low-me/smoothed-states/' 

smoothed_states = load_simulation()
smoothed_states_low_me = load_simulation(sim_dir=sim_dir_low_me)

smoothed_states['techdev'] = smoothed_states['Technology Level'] - smoothed_states['Technology Det']
smoothed_states_low_me['techdev'] = smoothed_states_low_me['Technology Level'] - smoothed_states_low_me['Technology Det']

mean = smoothed_states.groupby(level=0).mean()
q16 = smoothed_states.groupby(level=0).quantile(0.16)
q84 = smoothed_states.groupby(level=0).quantile(0.84)

mean_low_me = smoothed_states_low_me.groupby(level=0).mean()
q16_low_me = smoothed_states_low_me.groupby(level=0).quantile(0.16)
q84_low_me = smoothed_states_low_me.groupby(level=0).quantile(0.84)



cet_file = 'data/CET_AEJM_Fig5PanelC_Data.csv'
cet_data = p.read_csv(cet_file, index_col=0, parse_dates=True)
#cet_data.index = p.period_range(start='2008Q3', freq='Q', periods=cet_data.shape[0])
cet_data = cet_data.resample('Q',how=lambda x: x[-1])


#figure_defaults()
col = figure_defaults()

with saved_figure('tech_comp.pdf', nrows=1, ncols=2) as (fig, ax):


# Technology (bottom right)
    mean['techdev'].plot(ax=ax[0], linewidth=2)
    ax[0].fill_between(q16['techdev'].index, q16['techdev'], q84['techdev'], alpha=0.3)
    cet_data['CET Model'].plot(ax=ax[0],linestyle='dashed',color='purple')
    cet_data['Fernald (Util Adj.)'].plot(ax=ax[0],style='r:')
    ax[0].set_title(r'Technology Shock: Baseline ME',fontsize=16)
    despine()
#    ax[2,1].set_yticks(np.arange(-6,2,1))
#    ax[2,1].set_ylim([-6,1])
    ax[0].set_yticks(np.arange(-10,2,2))
    ax[0].set_ylim([-10,1])
##    ax[0].annotate('Fernald Utilization', xy=('2010Q1',0),fontsize=12,color='red')
##    ax[0].annotate('-adjusted TFP', xy=('2010Q2',-0.8),fontsize=12,color='red')
#   ax[0].annotate('CET TFP', xy=('2008Q2',-3.8),fontsize=12,color='green')
##    ax[0].annotate('CET TFP', xy=('2011Q3',-6),fontsize=12,color='purple')
    ax[0].grid(False)


    ax[0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0].tick_params('y',labelsize=14)
    ax[0].set_xlim(['2008','2014'])
    ax[0].set_xlabel(' ')


    mean_low_me['techdev'].plot(ax=ax[1], linewidth=2)
    ax[1].fill_between(q16_low_me['techdev'].index, q16_low_me['techdev'], q84_low_me['techdev'], alpha=0.3)
    cet_data['CET Model'].plot(ax=ax[1],linestyle='dashed',color='purple')
    cet_data['Fernald (Util Adj.)'].plot(ax=ax[1],style='r:')
    ax[1].set_title(r'Technology Shock: Low ME',fontsize=16)
    despine()

    ax[1].set_yticks(np.arange(-10,2,2))
    ax[1].set_ylim([-10,1])
##    ax[1].annotate('Fernald Utilization', xy=('2010Q1',0),fontsize=12,color='red')
##    ax[1].annotate('-adjusted TFP', xy=('2010Q2',-0.8),fontsize=12,color='red')
#   ax[1].annotate('CET TFP', xy=('2008Q2',-3.8),fontsize=12,color='green')
##    ax[1].annotate('CET TFP', xy=('2012Q2',-6),fontsize=12,color='purple')
    ax[1].grid(False)


    
    ax[1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1].tick_params('y',labelsize=14)
    ax[1].set_xlim(['2008','2014'])  
    ax[1].set_xlabel(' ')


    l1 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0])
    l2 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0],alpha=0.3)
    l3 = matplotlib.lines.Line2D([0],[0],linestyle='dashed',linewidth=2,color='purple')
    l4 = matplotlib.lines.Line2D([0],[0],linestyle='dotted',linewidth=2,color='red')
    
    fig.legend([l1,l2,l3,l4],['Model','68 Percent Credible Set','CET TFP','Fernald Utilization-adjusted TFP'],bbox_to_anchor=[0.5,0.05],loc='center',ncol=2,fontsize=10) 
    
    
    
    fig.tight_layout()
    #ax[2, 0].legend(['Smoothed_States', 'Only Liquidity Shocks', 'Only MEI Shocks', 'Only Tech. Shocks'], bbox_to_anchor=(0.5, 0.02), ncol=4)
    fig.subplots_adjust(bottom=0.15)
    #fig.subplots_adjust(hspace = .5, wspace = .4) # for 3x2 panels
    #fig.set_size_inches(5,7) # for 3x2 panels
    fig.subplots_adjust(hspace = .5, wspace = .3) 
    fig.set_size_inches(7, 3.5)

                        

