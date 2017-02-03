#### Used in aer2 figures

import numpy as np

from simulations_aer2 import load_simulation
from figures import saved_figure, despine, figure_defaults

import pandas as p
import matplotlib

from pandas import ExcelWriter

from matplotlib.ticker import MultipleLocator, FixedLocator


baseline = load_simulation()
nozlb = load_simulation(sim_dir='results/alt-sims/nozlb_2008')

gdp_diff = 100*(np.log(baseline.gdp) - np.log(nozlb.gdp))
gdp_mean = gdp_diff.groupby(level=0).mean()['2007Q4':]
gdp_q16 = gdp_diff.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
gdp_q84 = gdp_diff.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]
col = figure_defaults()

#nozlbq16 = nozlb.groupby(level=0).quantile(0.16)['Interest Rate']['2007Q4':]
#nozlbq84 = nozlb.groupby(level=0).quantile(0.84)['Interest Rate']['2007Q4':]

not_mean = baseline.groupby(level=0).mean()['Notional Rate']['2007Q4':]
not_q16 = baseline.groupby(level=0)['Notional Rate'].apply(lambda x: np.percentile(x, 16))['2007Q4':]
not_q84 = baseline.groupby(level=0)['Notional Rate'].apply(lambda x: np.percentile(x, 84))['2007Q4':]

inv_diff = 100*(np.log(baseline.inv) - np.log(nozlb.inv))
inv_mean = inv_diff.groupby(level=0).mean()['2007Q4':]
inv_q16 = inv_diff.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
inv_q84 = inv_diff.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]

cons_diff = 100*(np.log(baseline.cc) - np.log(nozlb.cc))
cons_mean = cons_diff.groupby(level=0).mean()['2007Q4':]
cons_q16 = cons_diff.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
cons_q84 = cons_diff.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]


######### Low Measurement Error
baseline_low_me = load_simulation_old(sim_dir=sim_dir_low_me)
nozlb_low_me = load_simulation_old('nozlb_2008_',sim_dir=sim_dir_low_me)

gdp_diff_low_me = 100*(np.log(baseline_low_me.gdp) - np.log(nozlb_low_me.gdp))
gdp_mean_low_me = gdp_diff_low_me.groupby(level=0).mean()['2007Q4':]
gdp_q16_low_me = gdp_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
gdp_q84_low_me = gdp_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]
col = figure_defaults()

#nozlbq16 = nozlb.groupby(level=0).quantile(0.16)['Interest Rate']['2007Q4':]
#nozlbq84 = nozlb.groupby(level=0).quantile(0.84)['Interest Rate']['2007Q4':]

not_mean_low_me = baseline_low_me.groupby(level=0).mean()['Notional Rate']['2007Q4':]
not_q16_low_me = baseline_low_me.groupby(level=0)['Notional Rate'].apply(lambda x: np.percentile(x, 16))['2007Q4':]
not_q84_low_me = baseline_low_me.groupby(level=0)['Notional Rate'].apply(lambda x: np.percentile(x, 84))['2007Q4':]

inv_diff_low_me = 100*(np.log(baseline_low_me.inv) - np.log(nozlb_low_me.inv))
inv_mean_low_me = inv_diff_low_me.groupby(level=0).mean()['2007Q4':]
inv_q16_low_me = inv_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
inv_q84_low_me = inv_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]

cons_diff_low_me = 100*(np.log(baseline_low_me.cc) - np.log(nozlb_low_me.cc))
cons_mean_low_me = cons_diff_low_me.groupby(level=0).mean()['2007Q4':]
cons_q16_low_me = cons_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 16))['2007Q4':]
cons_q84_low_me = cons_diff_low_me.groupby(level=0).apply(lambda x: np.percentile(x, 84))['2007Q4':]


with saved_figure('effect_of_zlb.pdf', nrows=2, ncols=2) as (fig, ax):
    
    not_mean.plot(ax=ax[0, 0], linewidth=2) 
    ax[0, 0].fill_between(not_q16.index, not_q16, not_q84, alpha=0.3)
    not_mean_low_me.plot(ax=ax[0, 0], linewidth=2, linestyle='dashed', color='red') 
    ax[0, 0].set_title(r'Notional Interest Rate',fontsize=16)
    ax[0, 0].set_xlim(['2007Q4','2014Q1'])
    ax[0, 0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0, 0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[0, 0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[0, 0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0, 0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[0, 0].tick_params('y',labelsize=14)
    ax[0, 0].grid(False)
    despine()
    
    gdp_mean.plot(ax=ax[0, 1], linewidth=2)
    ax[0, 1].fill_between(gdp_q16.index, gdp_q16, gdp_q84, alpha=0.3)
    gdp_mean_low_me.plot(ax=ax[0, 1], linewidth=2, linestyle='dashed', color='red') 
    ax[0, 1].set_title(r'Output',fontsize=16)
    ax[0, 1].set_xlim(['2007Q4','2014Q1'])
    ax[0, 1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0, 1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[0, 1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[0, 1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0, 1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[0, 1].tick_params('y',labelsize=14)
    ax[0, 1].set_yticklabels(np.arange(-4,2,1))
    ax[0 ,1].set_ylim([-4,1])
    ax[0, 1].grid(False)
    
    cons_mean.plot(ax=ax[1, 0], linewidth=2)
    ax[1, 0].fill_between(cons_q16.index, cons_q16, cons_q84, alpha=0.3)
    cons_mean_low_me.plot(ax=ax[1, 0], linewidth=2, linestyle='dashed', color='red') 
    ax[1, 0].set_title(r'Consumption',fontsize=16)
    ax[1, 0].set_xlim(['2007Q4','2014Q1'])
    ax[1, 0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1, 0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[1, 0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[1, 0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1, 0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[1, 0].tick_params('y',labelsize=14)
    ax[1, 0].grid(False)
    despine()
    
    inv_mean.plot(ax=ax[1, 1], linewidth=2)
    ax[1, 1].fill_between(inv_q16.index, inv_q16, inv_q84, alpha=0.3)
    inv_mean_low_me.plot(ax=ax[1, 1], linewidth=2, linestyle='dashed', color='red') 
    ax[1, 1].set_title(r'Investment',fontsize=16)
    despine()
    ax[1, 1].set_xlim(['2007Q4','2014Q1'])
    ax[1, 1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1, 1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[1, 1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[1, 1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1, 1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[1, 1].tick_params('y',labelsize=14)
    ax[1, 1].grid(False)

    #fig.tight_layout()
    fig.subplots_adjust(hspace=0.5,wspace=0.3)

    
    l1 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0])
    l2 = matplotlib.lines.Line2D([0],[0],linestyle='dashed',linewidth=2,color='red')
    l3 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0],alpha=0.3)
    
    fig.legend([l1,l2,l3],['Baseline ME', 'Low ME', '68 Percent Credible Set'],bbox_to_anchor=[0.43,0.02],loc='center',ncol=3,fontsize=10) 

    #ax.axhline(0, color='darkgray', linewidth=1)
    # xlabels = p.PeriodIndex(range(1985, 2016, 5), freq='A')
    # ax.set_xticklabels(xlabels)





# with saved_figure('effect_of_zlb_cont.pdf', nrows=2, ncols=2) as (fig, ax):
    

#     cons_mean.plot(ax=ax[0], linewidth=2)
#     ax[0].fill_between(cons_q16.index, cons_q16, cons_q84, alpha=0.3)
#     ax[0].set_title(r'Consumption')
#     despine()


#     fig.set_size_inches(12, 4)

