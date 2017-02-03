import numpy as np
import pandas as p

from model import varnames, shocknames, paranames

date_index = p.period_range(freq='Q', start='1982Q4', periods=126)

glss_names = ['beta','pibar','gz','psil','gamma','sigmal','phi','phiw',
              'ep','epw','ap','aw','bw','lamhp','alpha','delta','phii',
              'sigmaa','gam_rs','gam_dp','gamxhp','gamdy','shrgy',
              'sdevtech','rhog','sdevg','rhoinv','sdevinv','rholiq',
              'sdevliq','rhoint','sdevint','rhoelast','sdevelast',
              'rhoelastw','sdevelastw','meas_ygr', 'meas_cgr', 'meas_igr', 
              'meas_lab', 'meas_wgr', 'meas_infl', 'meas_nomr']



def from_matt(df):
    """Converts Matt's to my parameterization"""
    df['beta_tr'] = 100*(df['beta']**-1 - 1)
    df['pibar_tr'] = 100*(df['pibar'] - 1)
    df['gz_tr'] = 100*np.log(df['gz'])
  
    df['aw_tr'] =  1-df['aw']
    df['ap_tr'] =  1-df['ap']
    df['bw'] =  df['aw']
    df['sdevliq_tr'] =   df['sdevliq'] * 100.0
    df['sdevinv_tr'] =   df['sdevinv'] * 100.0
    df['sdevtech_tr'] =  df['sdevtech'] * 100.0
    df['sdevint_tr'] =  df['sdevint'] * 100.0
    df['sdevg_tr'] =  df['sdevg'] * 100.0
    df['sdevelast_tr'] =  df['sdevelast'] * 100.0
    df['sdevelastw_tr'] =  df['sdevelastw'] * 100.0

    return df


def to_matt(df):
    """Convert's my parameterization to Matt's"""
    df['gz'] = np.exp(df['gz_tr']/100.0)
    df['delta'] = 0.025
    df['shrgy'] = 0.2
    df['rholiq'] = 0.85
    df['ep_tr'] = 0.20
    df['epw_tr'] = 0.20
    df['rhoint'] = 0.0

    df['rhoelast'] =  0.0
    df['sdevelast_tr'] =   0.0
    df['rhoelastw'] =   0.0
    df['sdevelastw_tr'] =  0.0

    df['lamhp'] =  1600.00000000
    df['psil'] =  1.00000000
    df['beta'] =  1.0/(df['beta_tr']/100.0 + 1.0)
    df['pibar'] =  1.0/100*df['pibar_tr'] + 1
    df['ep'] =  (1.0+df['ep_tr'])/df['ep_tr']
    df['epw'] =  (1.0+df['epw_tr'])/df['epw_tr']
    df['aw'] =  1-df['aw_tr']
    df['ap'] =  1-df['ap_tr']
    df['bw'] =  df['aw']
    df['sdevliq'] =   df['sdevliq_tr'] / 100.0
    df['sdevinv'] =   df['sdevinv_tr'] / 100.0
    df['sdevtech'] =  df['sdevtech_tr'] / 100.0
    df['sdevint'] =  df['sdevint_tr'] / 100.0
    df['sdevg'] =  df['sdevg_tr'] / 100.0
    df['sdevelast'] =  df['sdevelast_tr'] / 100.0
    df['sdevelastw'] =  df['sdevelastw_tr'] / 100.0

    df['meas_ygr']= 0.00128825 * 2.5
    df['meas_cgr']= 0.00102312 * 2.5
    df['meas_igr']= 0.00485487 * 2.5
    df['meas_lab']= 0.0
    df['meas_wgr']= 0.0 
    df['meas_infl']= 0.0004877 * 2.5
    df['meas_nomr']= 0.0014285 * 2.5

    df['rholiq'] = 0.85

    return df


sim_dir = 'results/by-draw/'

#from tqdm import tqdm

import json

def load_simulation(date_index=date_index, sim_dir=sim_dir):
    def load_sim_i(i=0):

        ofile = sim_dir+'output{:04d}.json'.format(i)

        try:
            with open(ofile, 'r') as f:
                output = json.load(f)
        except:
            return 
     
        parasim =  p.DataFrame(output['input']['p0']).transpose()
        parasim.columns = glss_names
     
        parasim.gamtil = parasim.gamma/parasim.gz
        parasim.mc = (parasim.ep-1.0)/parasim.ep
        parasim.k2yrat = ((parasim.mc*parasim.alpha)/(parasim.gz/parasim.beta-(1.0-parasim.delta)))*parasim.gz
        parasim.shriy = (1.0-(1.0-parasim.delta)/parasim.gz)*parasim.k2yrat
        parasim.shrcy = 1.0-parasim.shrgy-parasim.shriy
        parasim.labss = ( ((parasim.epw-1.0)/parasim.epw)*(1.0-parasim.alpha)*(1.0-parasim.beta*parasim.gamtil)*((parasim.ep-1.0)/parasim.ep)*(1.0/(parasim.psil*(1.0-parasim.gamtil)))*(1.0/parasim.shrcy) )**(1.0/(parasim.sigmal+1.0))
     
        parasim['sig_uncon_liq'] = np.sqrt(parasim.sdevliq**2/(1-parasim.rholiq**2))
        parasim['sig_uncon_inv'] = np.sqrt(parasim.sdevinv**2/(1-parasim.rhoinv**2))
        parasim['sig_uncon_demand'] = np.sqrt(parasim.sdevg**2/(1-parasim.rhog**2))
        parasim['sig_uncon_mp'] = np.sqrt(parasim.sdevint**2)


        sims = [output['output']['sim_{:03d}'.format(i+1)]['smoothed_states'] 
                for i in range(10)]
        dfs = [p.DataFrame(sim, index=date_index).assign(sim=i) for i, sim in enumerate(sims)]
        df = p.concat(dfs)
        df = df.rename(columns=dict(zip(['endogvar_{:02d}'.format(i+1) for i in range(28)], 
                                        varnames[:28])))
        df['draw'] = i

        growth_rates = df.groupby(['sim']).apply(lambda x: 100*np.log(x).diff())

        df['Output Growth'] = growth_rates['gdp'] + 100*df.techshk + 100*np.log(parasim.gz.values)
        df['Consumption Growth'] = growth_rates['cc'] + 100*df.techshk + 100*np.log(parasim.gz.values)
        df['Investment Growth'] = growth_rates['inv'] + 100*df.techshk + 100*np.log(parasim.gz.values)
        df['Technology Growth'] = 100*df.techshk + 100*np.log(parasim.gz.values)
	#df['Technology Trend'] = 0.0

        df['Technology Trend'] = df['Technology Growth']
        df.loc[df.index > '2007Q4','Technology Trend'] = 100*np.log(parasim.gz.values)

        df['Inflation'] = 400*np.log(df['dp'])
        df['Interest Rate'] = 400*np.log(df['nomr'])
        df['Notional Rate'] = 400*np.log(df['notr'])
        df['Technology Level'] = df.groupby(['draw','sim']).cumsum()['Technology Growth']
        df['Technology Det'] = df.groupby(['draw','sim']).cumsum()['Technology Trend']
        
        df['xgap'] = parasim.alpha.values * np.log(df['util']).values + (1-parasim.alpha.values) * (np.log(df['lab'].values) -np.log(parasim.labss.values))
        df['scaled_liq'] = df['liqshk'] / parasim.sig_uncon_liq.values
        df['scaled_inv'] = df['invshk'] / parasim.sig_uncon_inv.values
        df['scaled_mp'] = df['intshk'] / parasim.sig_uncon_mp.values
        df['scaled_demand'] = df['gshk'] / parasim.sig_uncon_demand.values

        return df

    smoothed_states = p.concat([load_sim_i(i) for i in range(2000)])

    return smoothed_states
    
from model import varnames, shocknames, paranames

date_index = p.period_range(freq='Q', start='1982Q4', periods=126)

glss_names = ['beta','pibar','gz','psil','gamma','sigmal','phi','phiw',
              'ep','epw','ap','aw','bw','lamhp','alpha','delta','phii',
              'sigmaa','gam_rs','gam_dp','gamxhp','gamdy','shrgy',
              'sdevtech','rhog','sdevg','rhoinv','sdevinv','rholiq',
              'sdevliq','rhoint','sdevint','rhoelast','sdevelast',
              'rhoelastw','sdevelastw','meas_ygr', 'meas_cgr', 'meas_igr', 
              'meas_lab', 'meas_wgr', 'meas_infl', 'meas_nomr']


