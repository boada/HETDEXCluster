import h5py as hdf
from halo_handler import find_indices
from astLib import astStats
import numpy as np
from calc_cluster_props import findLOSV, findClusterCenterRedshift

f = hdf.File('out1204878.hdf5', 'r')
dset = f[f.keys()[2]]
data = dset.value

# now we need to make a mask for the data
mask = data['M200'] >= 1e14

# we'll use the mask to make all the changes and then consolidate back.
dataMasked = data[mask]
hids = np.unique(dataMasked['HALOID'])

halos = find_indices(dataMasked['HALOID'], hids)

for h in halos:
    if len(h) < 2:
        dataMasked['CLUSZ'][h] = np.mean(dataMasked['Z'][h])
    else:
        # find the cluster redshifts
        #dataMasked[h] = findClusterCenterRedshift(dataMasked['Z'][h])

        # find the LOSV
        #dataMasked[h] = findLOSV(dataMasked[h])

        # find the LOSVD
        dataMasked['LOSVD'][h] =\
        astStats.biweightScale_test(dataMasked['LOSV'][h], tuningConstant=9.0)

# update the data and file
#data[mask] = dataMasked
#dset[...] = data
#f.close()

