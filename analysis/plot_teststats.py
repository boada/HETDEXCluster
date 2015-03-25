import h5py as hdf
import pylab as pyl
from halo_handler import find_indices

f = hdf.File('./../data/halos/Aardvark_v1.0_halos_r1_rotated.4.hdf5')
#f = hdf.File('/home/boada/scratch/halos/Aardvark_v1.0_halos_r1_rotated.4.hdf5')
dset = f[f.keys()[0]]
halo = dset['HALOID', 'Z', 'VRMS', 'R200', 'M200']
f.close()

f = hdf.File('./stat_testing281041.hdf5', 'r')
dset = f[f.keys()[0]]
result = dset.value

# Only big clusters:
x = (halo['M200'] >= 1e13)
halo = halo[x]

halo = pyl.sort(halo, order = 'HALOID')
result = pyl.sort(result, axis=0)

inds = find_indices(halo['HALOID'], result[:,0])
halo = halo[inds]
inds = find_indices(result[:,0], halo['HALOID'])
result = result[inds]

pyl.scatter(result[:,1], (result[:,2] -
    halo['VRMS']/pyl.sqrt(3))/halo['VRMS']/pyl.sqrt(3))

pyl.show()


