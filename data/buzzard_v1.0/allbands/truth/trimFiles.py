import h5py as hdf
import numpy as np
from glob import glob

files = glob('*.hdf5')

for fname in files:
    print fname
    with hdf.File(fname, 'r') as f:
        key = f.keys()[0]
        dset = f[key]
        data = dset.value
        mask = data['Z'] < 0.5
    with hdf.File(fname.replace('.hdf5', '_trimmed.hdf5'), 'w') as f:
        f[key] = data[mask]

