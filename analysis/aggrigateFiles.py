import h5py as hdf
from glob import glob
import numpy as np

files = glob('*allGalaxies_DSresult.hdf5')

for f in files:
    with hdf.File(f, 'r') as hdfFile:
        dset = hdfFile[hdfFile.keys()[0]]
        data = dset.value

    print data.shape
    print data['HALOID'][:4]
    try:
        result = np.append(result, data)
    except NameError:
        result = data


#with hdf.File('out1204878_allGalaxies_DSresult.hdf5', 'w') as f:
#    f['ds_complete'] = result
#    f.flush()

