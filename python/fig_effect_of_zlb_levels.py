#### Used in aer2 figures
import numpy as np

from model import data
from simulations_aer2 import load_simulation
from figures import saved_figure, despine, figure_defaults

import pandas as p

from pandas import ExcelWriter

from matplotlib.ticker import MultipleLocator, FixedLocator

baseline = load_simulation()
nozlb = load_simulation(sim_dir='results/alt-sims/nozlb_2008')

init_cond = '2007Q4'
end_cond = '2014'
ss_inflation = 4*0.62 #400*np.log(parasim.mean().pibar)
ss_growth = 0.50 #100*np.log(parasim.mean().gz)

baseline['sspi'] = ss_inflation
baseline['ssy'] = ss_growth

#baseline.ix[init_cond]['Output Growth'] = 100.0
#nozlb.ix[init_cond]['Output Growth'] = 100.0

baseline.loc[init_cond,'Output Growth'] = 100.0
nozlb.loc[init_cond,'Output Growth'] = 100.0

mean_gdp = (baseline.groupby(level=0).mean()['Output Growth'][init_cond:end_cond]).cumsum()
mean_gdp_nozlb = (nozlb.groupby(level=0).mean()['Output Growth'][init_cond:end_cond]).cumsum()


data['Output Growth'][init_cond] = 100.0


baseline.loc[init_cond,'Inflation'] = 400.0
baseline.loc[init_cond,'sspi'] = 400.0
nozlb.loc[init_cond,'Inflation'] = 400.0

ss_inflation = (baseline.groupby(level=0).mean()['sspi'][init_cond:end_cond]/4.0).cumsum()
mean_pi = (baseline.groupby(level=0).mean()['Inflation'][init_cond:end_cond]/4.0).cumsum()
mean_pi_nozlb = (nozlb.groupby(level=0).mean()['Inflation'][init_cond:end_cond]/4.0).cumsum()


# 508 
# writer = ExcelWriter('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/effect_of_zlb_levels.xlsx')
# mean_gdp.to_frame(name='Output Growth').to_excel(writer,sheet_name='mean_gdp')
# mean_gdp_nozlb.to_frame(name='Output Growth').to_excel(writer,sheet_name='mean_gdp_nozlb')
# mean_pi.to_frame(name='Inflation').to_excel(writer,sheet_name='mean_pi')
# mean_pi_nozlb.to_frame(name='Inflation').to_excel(writer,sheet_name='mean_pi_nozlb')
# writer.save()



gdp = ((baseline.ix[init_cond:].groupby(level=1)).cumsum()['Output Growth'] - (nozlb.ix[init_cond:].groupby(level=1)).cumsum()['Output Growth']) / (baseline.ix[init_cond:].groupby(level=1)).cumsum()['Output Growth']

print((100*gdp).quantile([0.05, 0.16, 0.50, 0.84, 0.95]))
fjdkslfds

zlb_cost = p.DataFrame([mean_gdp - 100,  mean_gdp_nozlb - 100]).T
print "                            {:5s}     {:5s}".format('ZLB', 'NOZLB')
print "          2009Q2 Loss :   {:5.2f}     {:5.2f} ".format(*zlb_cost.ix['2009Q2'])
print "2008-2012 Avg Ann Loss:   {:5.2f}     {:5.2f} ".format(*zlb_cost.cumsum().ix['2012Q4']/(5*4))
fjdsklf
col = figure_defaults()

with saved_figure('effect_of_zlb_levels.pdf', ncols=2) as (fig, ax):
    
    #ss_growth.plot(ax=ax[0], linewidth=1)
    #ss_growth2.plot(ax=ax[0], linewidth=1)
    #data['Output Growth'][init_cond:].cumsum().plot(ax=ax[0])
    mean_gdp.plot(ax=ax[0])
    mean_gdp_nozlb.plot(ax=ax[0], linestyle='dashed')
    ax[0].set_title('Output',fontsize=16)
    ax[0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
    
    ax[0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])#,'2014Q3'])
    ax[0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])#,'2014'])
    
    ax[0].tick_params('y',labelsize=14)
    ax[0].set_xlim(['2008Q1','2014Q1'])#'2015Q1'])
    ax[0].annotate('Without ZLB', xy=('2008Q4',100),fontsize=14,color=col[1])
    ax[0].annotate('With ZLB', xy=('2010Q2',95),fontsize=14,color=col[0])
    
    ax[0].grid(False)

    despine()

    #ss_inflation.plot(ax=ax[1], linewidth=1)
    mean_pi.plot(ax=ax[1])
    mean_pi_nozlb.plot(ax=ax[1], linestyle='dashed')
    ax[1].set_title('Price Level',fontsize=16)
    ax[1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
    ax[1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])#,'2014Q3'])
    ax[1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])#,'2014'])
    ax[1].tick_params('y',labelsize=14)
    ax[1].set_xlim(['2008Q1','2014Q1'])#'2015Q1'])

    ax[1].grid(False)

    despine()

    fig.set_size_inches(8, 4)xcv
    fig.subplots_adjust(wspace = .25)
