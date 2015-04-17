import h5py as hdf
from halo_handler import find_indices
from astLib import astStats
f = hdf.File('out1204878.hdf5', 'r')
dset = f[f.keys()[2]]
data = dset.value
x = where(data['M200'] >= 1e14)
data = data[x[0]]
hid = unique(data['HALOID'])
halos = find_indices(data['HALOID'], hid)

for h in halos:
    #data['CLUSZ'][h] = astStats.biweightLocation(data['Z'][h],
    #        tuningConstant=6.0)
    #data['CLUSZ'][h] = mean(data['Z'][h])
