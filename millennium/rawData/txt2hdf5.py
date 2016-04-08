import numpy as np
import h5py as hdf
from glob import glob

halos = glob('halos*.txt')
galaxies = glob('galaxies*.txt')

for f in halos:
    d = np.genfromtxt(f, names=True, dtype=None)
    try:
        data = np.append(haloData, d)
    except NameError:
        haloData = d

for f in galaxies:
    d = np.genfromtxt(f, names=True, dtype=None)
    try:
        data = np.append(galData, d)
    except NameError:
        galData = d



with hdf.File('halos.hdf5', 'w') as f:
    f['halos'] = haloData

mask = galData['g_sdss'] == 99
with hdf.File('galaxies.hdf5', 'w') as f:
    f['galaxies'] = galData[~mask]

