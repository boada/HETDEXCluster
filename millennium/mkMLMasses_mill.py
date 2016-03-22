import numpy as np
import h5py as hdf
from sklearn.ensemble import RandomForestRegressor

with hdf.File('../analysis/result_targetedRealistic.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
                            'LOSVD_err', 'MASS']

    X = np.log10(target['M200c'])
    y = np.column_stack([np.log10(target['LOSVD']), target['ZSPEC'],
                            target['NGAL']])
