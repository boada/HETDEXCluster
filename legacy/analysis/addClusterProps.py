import h5py as hdf
from addHaloInfo import find_indices
from astLib import astStats
import numpy as np
from calc_cluster_props import findLOSV, findClusterCenterRedshift,\
    calc_mass_Evrard

def updateArray(data):
    from numpy.lib import recfunctions as rfns
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    print 'update array...'
    newData = -np.ones(len(data))
    data = rfns.append_fields(data, ['CLUSZ', 'LOSV', 'LOSVD', 'MASS'],
        [newData, newData, newData, newData], dtypes='>f4', usemask=False)

    return data

# silly stuff
import sys
import time
def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

spinner = spinning_cursor()

with hdf.File('out1204878_halo.hdf5', 'r') as f:
#with hdf.File('out1204878_allGalaxies.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

#data = updateArray(data)

# now we need to make a mask for the data -- HETDEX DEPTH!!!
#mask1 = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.5)
#mask2 = (data['g'] < 22.) | (data['Oii'] > 3.5)
#mask = mask1 & mask2

# we'll use the mask to make all the changes and then consolidate back.
#dataMasked = data[mask]
dataMasked = data

hids = np.unique(dataMasked['HALOID'])

halos = find_indices(dataMasked['HALOID'], hids)

for i,h in enumerate(halos):

    if len(h) < 5:
        dataMasked['CLUSZ'][h] = -1.
        dataMasked['LOSV'][h] = -1.
        dataMasked['LOSVD'][h] = -1.

    elif 5 <= len(h) < 15:
        dataMasked[h] = findClusterCenterRedshift(dataMasked[h])
        dataMasked[h] = findLOSV(dataMasked[h])
        dataMasked['LOSVD'][h] =\
            astStats.gapperEstimator(dataMasked['LOSV'][h])

    elif 15 <= len(h):
        #find the cluster redshifts
        dataMasked[h] = findClusterCenterRedshift(dataMasked[h])
        #find the LOSV
        dataMasked[h] = findLOSV(dataMasked[h])

        # find the LOSVD
        dataMasked['LOSVD'][h] =\
            astStats.biweightScale_test(dataMasked['LOSV'][h],
                    tuningConstant=9.0)

    if not dataMasked['CLUSZ'][h][0] == -1.0:
        # finally the mass
        dataMasked['MASS'][h] = calc_mass_Evrard(dataMasked[h])
    else:
        pass

    sys.stdout.write(spinner.next())
    sys.stdout.flush()
    sys.stdout.write('\b')

# now we make another new file
with hdf.File('out1204878_complete.hdf5', 'w') as f:
#with hdf.File('out1204878_hetdex.hdf5', 'w') as f:
#with hdf.File('out1204878_allGalaxies_props.hdf5', 'w') as f:
#    data[mask] = dataMasked
    f['dset_complete'] = dataMasked
    f.flush()

