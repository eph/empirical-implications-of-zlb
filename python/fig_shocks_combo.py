#### Used in aer2 figures

from simulations_aer2 import load_simulation, date_index, sim_dir
from figures import saved_figure, despine, figure_defaults

import numpy as np
import pandas as p

from pandas import ExcelWriter

from matplotlib.ticker import MultipleLocator, FixedLocator


baseline = load_simulation()


mean = baseline.groupby(level=0).mean()['2008':]
q16 = baseline.groupby(level=0).quantile(0.16)['2008':]
q84 = baseline.groupby(level=0).quantile(0.84)['2008':]

#shocks = ['scaled_liq', 'scaled_inv', 'techdev']

#names = {'scaled_liq': 'Risk Premium Shock', 
#         'scaled_inv': 'MEI Shock', 
#         'techdev': 'Technology Shock'}


ebp = p.read_csv('data/GLZ_variables.csv', index_col=0, parse_dates=True).ebp_oa
ebpq = ebp.resample('Q', how=lambda x: x[-1])
ebpq = ebpq / ebpq.std()

baa = p.read_csv('data/Baaspread10y.csv', index_col = 0, parse_dates=True).baa
baaq = baa.resample('Q',how=lambda x: x[-1])

erp = [3.7798,  3.6601,   3.9595,   3.9612,   3.8119,   3.7339,   3.8315,   4.9372,   4.8648,   4.7810,   5.8417,   5.5786,   7.4788,   7.7432,   8.3548,   8.0830,   7.4911,   7.8343,   8.6081,   8.4127,   9.4940,  11.5793,  12.8545,  12.9051,  12.5041,  13.0734,  13.6974,  12.1652,  11.6023,  11.3401,  10.6794,  10.1015,   9.8318,   9.2482,   9.1041,   8.8669,   8.5851,   8.7962,   8.2809,   8.0348,   8.5987,   8.9522,   9.1013,   9.1247,   8.8152,   8.6266,   8.3927,   8.0683,   8.0297,   7.7565,   7.9376,   8.1317,   8.1017,   8.4215,   8.4493,   9.4906,   9.6438,   9.7760,   9.6719,   9.5697,   9.4257,   9.0696,   8.8413,   9.1381,   9.4851,   9.5763,   9.7325,   9.4708,   9.2138,   9.6695,   9.9707,   9.8267,   9.6674,   9.4558,   9.2933,   9.5280,   9.1225,   9.2017,   9.1923,   9.1795,   9.1234,   9.0744,   8.7467,   8.6200,   8.9146,   8.8910,   8.6727,   8.9890,   8.8787,   8.6189,   8.7492,   8.8183,   8.6794,   9.1658,   8.7191,   8.6609,   9.1134,   8.8676,   8.8663,   8.9956,   8.9457,   8.9604,   9.0846,   9.2850,   9.8057,   9.6808,   9.3658]

erp = p.Series(erp, index=p.period_range(freq='M', start='2007-1', periods=len(erp))).resample('Q')

eerp = p.read_csv('data/EquityRiskPremiumTR.csv', index_col = 0, parse_dates=True).eerp
eerp = eerp['2007':'2014']
eerpq = eerp.resample('Q',how='mean')


#rp_innovation_mean = p.read_csv('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/rp_innovation.csv',index_col = 0, header=None, names = ['x','rp'])

# RP Innovation
from os import listdir; 
output_dir = sim_dir

output_list = [f for f in listdir(output_dir) 
               if f.startswith('output')]

import json


#filtered_states = p.DataFrame()
#smoothed_states = p.DataFrame()
smoothed_shocks = p.DataFrame()
#filtered_shocks = p.DataFrame()
j = 0
from tqdm import tqdm 
for output in tqdm(output_list):
    
    data = json.load(open(output_dir+'/'+output))

    for sim in data['output']:

        # smooth_i = p.DataFrame(data['output'][sim]['smoothed_states'])
        # filter_i = p.DataFrame(data['output'][sim]['filtered_states'])

        # smooth_i['trial'] = j
        # filter_i['trial'] = j

        # smoothed_states = smoothed_states.append(smooth_i)
        # filtered_states = filtered_states.append(filter_i)

        smooth_i = p.DataFrame(data['output'][sim]['smoothed_shocks'])
#        filter_i = p.DataFrame(data['output'][sim]['filtered_shocks'])

        smooth_i['trial'] = j
#        filter_i['trial'] = j

        smoothed_shocks = smoothed_shocks.append(smooth_i)


        j = j+1

#print(output_list)

date_index = p.period_range(freq='Q', start='1982Q4', periods=126)
shocks_mean = smoothed_shocks.groupby(level=0).mean()
shocks_q16 = smoothed_shocks.groupby(level=0).quantile(0.16)
shocks_q84 = smoothed_shocks.groupby(level=0).quantile(0.84)

rp_mean = p.Series(shocks_mean['exogvar_01'])
rp_mean.index = date_index

rp_q16 = p.Series(shocks_q16['exogvar_01'])
rp_q16.index = date_index

rp_q84 = p.Series(shocks_q84['exogvar_01'])
rp_q84.index = date_index

# Equity premium
premium_dir='results/alt-sims/'

premium_list = [f for f in listdir(premium_dir) 
               if f.startswith('equity_premiumoutput')]

premium = p.DataFrame()
j = 0

for output_p in tqdm(premium_list):
    
    data = json.load(open(premium_dir+'/'+output_p))

    for sim in data['output']:

        premium_i = p.DataFrame(data['output'][sim]['smoothed_states'])
        #premium_i['trial'] = j

        premium = premium.append(premium_i)

        j = j+1

premium = 400*np.log(premium)

premium_mean = premium.groupby(level=0).mean()
premium_q16 = premium.groupby(level=0).quantile(0.16)
premium_q84 = premium.groupby(level=0).quantile(0.84)

ep_mean = p.Series(premium_mean['endogvar_29'])
ep_mean.index = date_index

ep_q16 = p.Series(premium_q16['endogvar_29'])
ep_q16.index = date_index

ep_q84 = p.Series(premium_q84['endogvar_29'])
ep_q84.index = date_index




figure_defaults()


with saved_figure('shocks_combo.pdf', nrows=2, ncols=2) as (fig, ax):
    
    # Risk Premium Shock (top left)
    mean['scaled_liq'].plot(ax=ax[0, 0], linewidth=2)
    ax[0, 0].fill_between(q16['scaled_liq'].index, q16['scaled_liq'], q84['scaled_liq'], alpha=0.3)
    ebpq['2008':].plot(ax=ax[0, 0], linestyle='dashed', color='green')
    baaq['2008':].plot(ax=ax[0, 0],linestyle='dotted',color='red')
    ax[0, 0].set_title(r'Risk Premium Shock',fontsize=16)
    ax[0, 0].set_xlim(['2008Q1','2014Q1'])
    ax[0, 0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0, 0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[0, 0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[0, 0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0, 0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[0, 0].tick_params('y',labelsize=14)
    ax[0, 0].set_yticklabels(np.arange(-2,7,1))
    ax[0 ,0].set_ylim([-2,6])
    ax[0, 0].annotate('GZ Spread', xy=('2009Q1',-1.65),color='green',fontsize=12)
    ax[0, 0].annotate('BAA Corporate Bond Spread', xy=('2009Q3',4),color='red',fontsize=12)
    ax[0, 0].set_xlabel('')
    ax[0, 0].grid(False)
#    ax[0, 0].spines['left'].set_color('k')
#    ax[0, 0].spines['bottom'].set_color('k')



    # MEI (top right)
    mean['scaled_inv'].plot(ax=ax[0, 1], linewidth=2)
    ax[0, 1].fill_between(q16['scaled_inv'].index, q16['scaled_inv'], q84['scaled_inv'], alpha=0.3)
    ax[0, 1].set_title(r'MEI Shock',fontsize=16)
    ax[0, 1].set_xlim(['2008Q1','2014Q1'])
    ax[0, 1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[0, 1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[0, 1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[0, 1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[0, 1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[0, 1].tick_params('y',labelsize=14)
    ax[0, 1].set_yticks(np.arange(-3,3,1))
    ax[0, 1].set_ylim([-3,2])
    ax[0, 1].grid(False)  
    
    despine()

    # Risk Premium Innovation (bottom left)
    rp_mean.plot(ax=ax[1,0],  linewidth=2)
    ax[1,0].fill_between(rp_q16.index, rp_q16, rp_q84, alpha=0.3)
    #mean['smoothed_shock_exogvar_01'].plot(ax=ax[1, 0], linewidth=2)
    ax[1,0].set_ylim([-1,5])
    ax[1,0].grid(False)
    ax[1,0].tick_params('y',labelsize=14)
    ax[1,0].set_title(r'Risk Premium Innovation',fontsize=16)
    ax[1,0].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1,0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[1,0].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[1,0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1,0].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[1,0].set_xlim(['2008Q1','2014Q1'])
    
    #rp_innovation_mean['rp'].plot(ax=ax[1,1])

    # Equity Premium (bottom right)
    #mean.premium.plot(ax=ax[0, 1], linewidth=2)
    ep_mean.plot(ax=ax[1, 1], linewidth=2)
    #ax[0, 1].fill_between(q16.premium.index, q16.premium, q84.premium, alpha=0.3)
    ax[1, 1].fill_between(ep_q16.index, ep_q16, ep_q84, alpha=0.3)
    erp['2008':'2014Q1'].plot(ax=ax[1, 1], linestyle='dashed', color='green')
    #eerpq['2008':].plot(ax=ax[0, 1], linestyle='dotted', color='red')
    ax[1, 1].set_title('Equity Premium',fontsize=16)
    ax[1, 1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[1, 1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=14,pad=5)
    ax[1, 1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[1, 1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=12,pad=5)
    ax[1, 1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[1, 1].tick_params('y',labelsize=14)
    ax[1, 1].set_xlim(['2008Q1','2014Q1'])
    ax[1, 1].annotate('Adrian et al. (2012)', xy=('2009Q4',11.5),color='green',fontsize=12)
    ax[1, 1].grid(False)
    ax[1, 1].set_ylim(0, 15)
    #ax[0, 1].annotate('Estimated Equity Risk Premium', xy=('2008Q3',3),color='red',fontsize=14)
    despine()
    fig.subplots_adjust(hspace=0.5,wspace=0.3)

