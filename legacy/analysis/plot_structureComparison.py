import h5py as hdf
import pylab as pyl
from addHaloInfo import find_indices

with hdf.File('out1204878_allGalaxies_DS.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    ds = dset.value
with hdf.File('out1204878_allGalaxies_props.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value
with hdf.File('out1204878_allGalaxies_ADSW.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    adsw = dset.value

# sort the data -- TODO:this should be already done in the files
mask = pyl.argsort(ds['HALOID'])
ds = ds[mask]
mask = pyl.argsort(adsw['HALOID'])
adsw = adsw[mask]

#hids = pyl.unique(data['HALOID'])
halos = pyl.array(find_indices(data['HALOID'], ds['HALOID']))

# masses
mass = pyl.array([data['MASS'][h][0] for h in halos])

# filter out bad values
mask = mass > 0
mass = mass[mask]
ds = ds[mask]
adsw = adsw[mask]


# now we are making the masks for the histograms
mask = adsw['A2'] > adsw['CRIT'][:,2]
mask = ds['DS'] < 0.01 && ds['DS'] >= 0


