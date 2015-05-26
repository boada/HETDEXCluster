import numpy as np
import h5py as hdf
from addHaloInfo import find_indices
from multiprocessing import Pool, cpu_count, current_process
from itertools import izip, repeat
from os import environ
import sys
from scipy.stats import anderson, shapiro

def mp_worker((h, data)):
    if len(h) < 10:
        return -1., -1, -1
    else:
        LOSV = data['LOSV'][h]
        A2, critical, _ = anderson(LOSV)
        w, p = shapiro(LOSV)

    return A2, critical, w

def mp_worker_wrapper(*args):
    return mp_worker(*args)

def start_process():
    print 'Starting', current_process().name

def main(start=None, end=None):
    with hdf.File('out1204878_allGalaxies_props.hdf5', 'r') as f:
        dset = f[f.keys()[0]]
        data = dset['LOSV', 'HALOID']

    hids = np.unique(data['HALOID'])
    halos = np.array(find_indices(data['HALOID'], hids))

    p = Pool(cpu_count(), maxtasksperchild=1, initializer=start_process)

    result = p.map(mp_worker_wrapper, izip(halos, repeat(data)),
            chunksize=500)

    p.close()
    p.join()

    finalResult = np.ones(hids.shape, dtype=[('HALOID', '>i8'),
        ('A2', '>f4'),  ('CRIT', '>(5,)f4'),  ('SW-P', '>f4')])
    finalResult['HALOID'] = hids
    for i, (A2, C, W) in enumerate(result):
        finalResult['A2'][i] = A2
        finalResult['CRIT'][i] = C
        finalResult['SW-P'][i] = W

    with hdf.File('out1204878_allGalaxies_ADSW.hdf5', 'w') as f:
        f['DSresult'] = finalResult
        f.flush()

    return finalResult
if __name__ == "__main__":
    main()
