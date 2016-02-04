import numpy as np
import h5py as hdf
from glob import glob

files = glob('*Realistic*')

for f in files:
    with hdf.File(f, 'r') as f2:
        dset = f2[f2.keys()[0]]
        print dset.file
        result_part = dset.value
        print result_part.size
    try:
        result = np.append(result, result_part)
    except NameError:
        result = result_part

with hdf.File('result_targetedRealistic.hdf5', 'w') as f:
    f['result_targetedReastlic'] = result

