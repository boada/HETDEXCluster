from multiprocessing import Pool
import h5py as hdf
import numpy as np
from halo_handler import find_indices
from calc_cluster_props import (updateArray, findClusterRedshift, findLOSV,
                                findLOSVDmcmc, calc_mass_Evrard)
import os


class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool(maxtasksperchild=10)

    def call(self, *args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()


def worker(pos, data):
    print("PID: %d \t Value: %d" % (os.getpid(), pos))
    data = updateArray(data)
    data = findClusterRedshift(data)
    #data = findSeperationSpatial(data, center)
    data = findLOSV(data)
    data, sigma_dist = findLOSVDmcmc(data)
    data = calc_mass_Evrard(data, A1D=1177, alpha=0.364)
    return pos, data, sigma_dist


def cb_func(result):
    (pos, data, sigma_dist) = result
    if pos % 1000 == 0:
        print(pos)
    results['IDX'][pos] = pos
    results['CLUSZ'][pos] = data['CLUSZ'][0]
    results['LOSVD'][pos] = data['LOSVD'][0]
    results['MASS'][pos] = data['MASS'][0]
    results['LOSVD_err'][pos] = data['LOSVD_err'][0]
    results['LOSVD_dist'][pos] = sigma_dist[:, 0]


if __name__ == "__main__":
    async_worker = AsyncFactory(worker, cb_func)

    with hdf.File('./halos.hdf5', 'r') as f:
        dset = f[list(f.keys())[0]]
        halo = dset.value
    with hdf.File('./galaxies_Oii.hdf5', 'r') as f:
        dset = f[list(f.keys())[0]]
        truth = dset['fofid', 'velx', 'g', 'redshift', 'galaxyid', 'Oii', 'Z']

    # this is the part that makes it realistic or not
    gmask = truth['g'] < 22
    Oiimask = truth['Oii'] > 3.5
    mask = gmask | Oiimask
    truth = truth[mask]

    hids, uniqueIdx = np.unique(halo['fofId'], return_index=True)

    # now find all the corresponding galaxies
    gals = find_indices(truth['fofid'], hids)

    # make the results container
    x = [i for i, g in enumerate(gals) if g.size >= 5]
    # make the results container
    results = np.zeros(
        (len(x), ),
        dtype=[('IDX', '>i4'), ('HALOID', '>i8'), ('ZSPEC', '>f4'),
               ('M200c', '>f4'), ('CLUSZ', '>f4'), ('LOSVD', '>f4'),
               ('MASS', '>f4'), ('NGAL', '>i4'), ('LOSVD_err', '>f4',
                                                  (2, )), ('LOSVD_dist', '>f4',
                                                           (10000, ))])

    print(('do work', len(x), 'clusters to go!'))
    keepBad = False
    for j, i in enumerate(x):
        if gals[i].size >= 5:
            async_worker.call(j, truth[gals[i]])
            # update results array
            results['HALOID'][j] = halo['fofId'][uniqueIdx[i]]
            results['NGAL'][j] = gals[i].size
            results['ZSPEC'][j] = halo['redshift'][uniqueIdx[i]]
            results['M200c'][j] = halo['m_crit200'][uniqueIdx[i]] * 1e10 / 0.73

    async_worker.wait()

    try:
        #os.remove('result_targetedPerfect.hdf5')
        os.remove('result_targetedRealistic.hdf5')
    except OSError:
        pass
    #with hdf.File('result_targetedPerfect.hdf5', 'w') as f:
    #    f['result_targetedPerfect'] = results
    with hdf.File('result_targetedRealistic.hdf5', 'w') as f:
        f['result_targetedRealistic'] = results
