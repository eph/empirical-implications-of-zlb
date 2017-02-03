from model import data
from simulations_aer2 import load_simulation
from figures import saved_figure, despine, figure_defaults

import pandas as p
import numpy as np
import matplotlib

from matplotlib.ticker import MultipleLocator, FixedLocator

from pandas import ExcelWriter

smoothed_states = load_simulation()

mean = smoothed_states.groupby(level=0).mean()
q16 = smoothed_states.groupby(level=0).quantile(0.16)
q84 = smoothed_states.groupby(level=0).quantile(0.84)

to_plot = ['Output Growth', 
           'Investment Growth', 
           'Consumption Growth', 
           'Inflation', 
           'Interest Rate', 
           'Notional Rate']

# 508
# writer = ExcelWriter('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/observable_fit.xlsx')
# data.to_excel(writer,sheet_name='data')
# mean.to_excel(writer,sheet_name='mean')
# q16.to_excel(writer,sheet_name='q16')
# q84.to_excel(writer,sheet_name='q84')
# writer.save()


col = figure_defaults()


#with saved_figure('/msu/res2/Shared_Projects/GLS/G_LS_S/latex_aer_revision/figures_aer2/observable_fit.pdf', nrows=3, ncols=2) as (fig, ax):
with saved_figure('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/figures_aer2/observable_fit.pdf', nrows=3, ncols=2) as (fig, ax):
    
    xlabels = p.PeriodIndex(range(1985, 2016, 5), freq='A')

    for axi, obs in zip(ax.reshape(-1), to_plot):

        if not(obs is 'Notional Rate'):
            data[obs].plot(ax=axi, linestyle='dashed', linewidth=2)
            mean[obs].plot(ax=axi, linewidth=2)
            axi.fill_between(q16[obs].index, q16[obs], q84[obs], alpha=0.3)
            axi.set_title(obs,fontsize=14)
            #axi.set_xticklabels(xlabels)

     
        else:

            mean[obs]['2008':].plot(ax=axi, linewidth=2)
            axi.fill_between(q16[obs]['2008':].index, q16[obs]['2008':], 
                             q84[obs]['2008':], alpha=0.3)
            axi.set_title(obs,fontsize=14) 
            
    
    for i in range(0,2):
        for j in range(0,2):
            ax[i,j].get_xaxis().set_minor_locator(MultipleLocator(20))
            ax[i,j].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
            ax[i,j].set_xticks(['1985','1990','1995','2000','2005','2010','2015'])
            ax[i,j].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=10,pad=5)
            ax[i,j].set_xticklabels(['1985','1990','1995','2000','2005','2010','2015'])
            ax[i,j].tick_params('y',labelsize=14)
            ax[i,j].set_xlim(['1985','2015'])  
            
    ax[2,0].get_xaxis().set_minor_locator(MultipleLocator(20))
    ax[2,0].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
    ax[2,0].set_xticks(['1985','1990','1995','2000','2005','2010','2015'])
    ax[2,0].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=10,pad=5)
    ax[2,0].set_xticklabels(['1985','1990','1995','2000','2005','2010','2015'])
    ax[2,0].tick_params('y',labelsize=14)
    ax[2,0].set_xlim(['1985','2015'])          
    
    
    ax[2,1].get_xaxis().set_minor_locator(MultipleLocator(4))
    ax[2,1].tick_params('x',length=5,width=1,which='minor',direction='in',labelsize=10,pad=5)
   # ax[2,1].set_xticklabels(p.PeriodIndex(range(2008,2015,1),freq='A'))
    ax[2,1].set_xticks(['2008Q3','2009Q3','2010Q3','2011Q3','2012Q3','2013Q3'])
    ax[2,1].set_xticklabels(['2008','2009','2010','2011','2012','2013'])
    ax[2,1].tick_params('x', length=0,width=0,which='major',direction='in',labelsize=10,pad=5)
    ax[2,1].tick_params('y',labelsize=14)
    ax[2,1].set_xlim(['2008','2014'])  
        
    ax[0,0].set_yticks(np.arange(-4,5,2))
    ax[0,0].set_ylim([-4,4])
    ax[0,0].grid(False)
        
    ax[0,1].set_yticks(np.arange(-15,11,5))
    ax[0,1].set_ylim([-15,10])
    ax[0,1].grid(False)        

    ax[1,0].set_yticks(np.arange(-3,3,1))
    ax[1,0].set_ylim([-3,2])
    ax[1,0].grid(False)
        
    ax[1,1].set_yticks(np.arange(-2,7,2))
    ax[1,1].set_ylim([-2,6])
    ax[1,1].grid(False)
        
    ax[2,0].set_yticks(np.arange(0,13,2))
    ax[2,0].set_ylim([0,12])
    ax[2,0].grid(False)
        
    ax[2,1].set_yticks(np.arange(-9,7,3))
    ax[2,1].set_ylim([-9,6])
    ax[2,1].grid(False)    
    
     
    despine()     
    
    #fig.legend(,to_plot,loc=(.5,.5),ncol=1)
    l1 = matplotlib.lines.Line2D([0],[0],linestyle='dashed',linewidth=2,color=col[0])
    l2 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[1])
    l3 = matplotlib.lines.Line2D([0],[0],linewidth=2,color=col[0],alpha=0.3)
    
    fig.legend([l1,l2,l3],['Data','Model','68 Percent Credible Set'],bbox_to_anchor=[0.5,0.02],loc='center',ncol=3,fontsize=13) 
        
    
    fig.tight_layout()
    fig.subplots_adjust(hspace = .5, wspace = .4)
    fig.set_size_inches(6,8)


