import h5py as hdf
from updateResults import find_indices
from astLib import astStats
import numpy as np
from calc_cluster_props import findLOSV, findClusterCenterRedshift, calc_mass_Saro

# silly stuff
import sys
import time
def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

spinner = spinning_cursor()

with hdf.File('out1204878_halo.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

    # now we need to make a mask for the data
    #mask = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.2)

    # we'll use the mask to make all the changes and then consolidate back.
    #dataMasked = data[mask]
    dataMasked = data
    hids = np.unique(dataMasked['HALOID'])

    halos = find_indices(dataMasked['HALOID'], hids)

    for i,h in enumerate(halos):
        if len(h) < 4:
            dataMasked['CLUSZ'][h] = np.mean(dataMasked['Z'][h])
            #find the LOSV
            dataMasked[h] = findLOSV(dataMasked[h])

            dataMasked['LOSVD'][h] = np.std(dataMasked['LOSV'][h], ddof=1.5)

        else:
            try:
                #find the cluster redshifts
                dataMasked[h] = findClusterCenterRedshift(dataMasked[h])
                #find the LOSV
                dataMasked[h] = findLOSV(dataMasked[h])

                # find the LOSVD
                dataMasked['LOSVD'][h] =\
                astStats.biweightScale_test(dataMasked['LOSV'][h],
                        tuningConstant=9.0)

            except ZeroDivisionError:
                badData = -np.ones(len(h))
                dataMasked['CLUSZ'][h] = badData
                dataMasked['LOSV'][h] = badData
                dataMasked['LOSVD'][h] = badData

        if not dataMasked['CLUSZ'][h][0] == -1.0:
            # finally the mass
            dataMasked['MASS'][h] = calc_mass_Saro(dataMasked[h])
        else:
            pass

        sys.stdout.write(spinner.next())
        sys.stdout.flush()
        sys.stdout.write('\b')

# now we make another new file
with hdf.File('out1204878_complete.hdf5', 'w') as f:
    #data[mask] = dataMasked
    f['dset_complete'] = dataMasked
    f.flush()

