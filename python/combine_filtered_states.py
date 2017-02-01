from __future__ import division

import sys

import pandas as p
import numpy as np

sys.path.append('/msu/scratch3/m1cjg01/aer_revision_ed/python/notebooks/') 

from model import varnames 

import matplotlib.pyplot as plt

nprocs = 200
npart = 150000
weights_file = 'test-output/t_{:03d}_rank_{:03d}_weights.txt'.format
states_file = 'test-output/t_{:03d}_rank_{:03d}_states.txt'.format
rs = []

from progressbar import ProgressBar
pb = ProgressBar()
for t in pb(range(1, 126)):
    arrays = [np.loadtxt(weights_file(t, i)) for i in range(nprocs)]
    weights = np.concatenate(arrays)
    assert len(weights) == npart
    npart_per_proc =  npart / nprocs 

    randi = np.random.choice(np.arange(npart), p=weights)

    proci, randi = divmod(randi, npart_per_proc)

    states = np.loadtxt(states_file(t, int(proci)))[i, :]

    rs.append(dict(zip(varnames, states)))

res = p.DataFrame(rs)