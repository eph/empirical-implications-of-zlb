import pandas as p
from simulations_aer2 import from_matt, glss_names

sim_dir =  'results/mcmc/'

paraest = p.DataFrame()
for i in range(5):
    parasim = p.read_csv(sim_dir + 'output{:i}_para_tmp.txt'.format(i), delim_whitespace=True, 
                         names=glss_names, index_col=False)
    parasim = parasim.ix[int(0.2*parasim.shape[0]):, :]


for i in range(2000):
    np.savetxt('results/thinned_posterior/parasim{:04d}.txt'.format(i), parasim.ix[::80].ix[i])

