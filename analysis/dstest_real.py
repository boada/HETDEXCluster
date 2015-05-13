import numpy as np
from sklearn.neighbors import NearestNeighbors
from calc_cluster_props import *
from astLib import astStats
from addHaloInfo import find_indices
from scipy import stats
from dstest import DStest


#with hdf.File('out1204878_hetdex.hdf5', 'r') as f:
with hdf.File('out1204878_complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

mask = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.5)
data = data[mask]

hids = np.unique(data['HALOID'])
halos = np.array(find_indices(data['HALOID'], hids))

blah = []
for h in halos[:10]:
    if len(h) < 10:
        pass
    else:
        cluster = np.column_stack([data['RA'][h], data['DEC'][h],
                data['Z'][h]])
        s = DStest(cluster, data['LOSV'][h], data['LOSVD'][h][0], shuffles=1000)

        blah.append([len(h), s])
#    print delta/len(cluster)
