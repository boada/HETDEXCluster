import pandas as pd
import h5py as hdf
from glob import glob

files = glob('*.hdf5')

for f in files:
    hdf5File = hdf.File(f, 'r')
    dset = hdf5File[hdf5File.keys()[0]]

    if 'df' not in vars():
        df = pd.DataFrame({'HALOID': dset['HALOID'], 'R200': dset['R200'], 'Z':
            dset['Z'], 'NGALS': dset['NGALS']})
    else:
        df = pd.concat([df, pd.DataFrame({'HALOID': dset['HALOID'], 'R200':
                            dset['R200'], 'Z': dset['Z'], 'NGALS':
                            dset['NGALS']})], ignore_index=True)

    hdf5File.close()
