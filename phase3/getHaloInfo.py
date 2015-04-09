import h5py as hdf
import numpy as np
from halo_handler import *

# load the result file
f = hdf.File('out1204878.hdf5', 'r')
dset = f[f.keys()[0]]

# load the coordinate information
RAmax = dset['RA'].max()
RAmin = dset['RA'].min()
DECmax = dset['DEC'].max()
DECmin = dset['DEC'].min()

# Get the unique haloids
haloids = np.unique(dset['HALOID'])

# find the tiles that overlap with our area
haloTiles = findHaloTile(RAmax, DECmax, RAmin, DECmin)

# loads the halo files

haloCat = mk_haloCatalog(haloTiles)



