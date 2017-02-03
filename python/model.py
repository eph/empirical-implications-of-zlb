import numpy as np

from dsge.DSGE import DSGE
import pandas as p 

mod_dir = 'python/'
mod_file = mod_dir + 'model.yaml'
mod = DSGE.read(mod_file)

GHLS = mod.compile_model()
GHLS.yy.index = p.period_range(freq='Q', start='1983Q1', periods=125)

texnames = mod['__data__']['tex']

prior_exo = [('rhog', 'sdevg_tr'), ('rhoinv', 'sdevinv_tr'), 
             ('rholiq', 'sdevliq_tr'), ('sdevtech_tr', 'sdevint_tr')]

prior_ss = [('rbar', 'pibar_tr'), ('gz_tr', 'alpha')]

prior_rule = [('gam_rs', 'gam_dp'), ('gamdy', 'gamxhp')]

prior_endo = [('gamma', 'sigmal'), ('sigmaa', 'phii'), 
              ('phi', 'ap_tr'), ('phiw', 'aw_tr')]


varnames = map(str, mod.variables)
shocknames = map(str, mod.shocks)
paranames = map(str, mod.parameters)


data = GHLS.yy

data['Output Growth'] = 100*data.ygr
data['Investment Growth'] = 100*data.igr
data['Consumption Growth'] = 100*data.cgr
data['Inflation'] = 400*data.infl 
data['Interest Rate'] = 400*data.nomr
data['Notional Rate'] = np.nan