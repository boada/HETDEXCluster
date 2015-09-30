from multiprocessing import Pool
import h5py as hdf
import numpy as np
from calc_cluster_props import *
import os

class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool()

    def call(self,*args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()

def worker(pos, data, center):
    #print "PID: %d \t Value: %d" % (os.getpid(), pos)
    data = updateArray(data)
    data = findClusterRedshift(data)
    #data = findSeperationSpatial(data, center)
    data = findLOSV(data)
    data = findLOSVD(data)
    data = findLOSVDgmm(data)
    data = calc_mass_Evrard(data, A1D = 1177, alpha = 0.364)
    return pos, data

def cb_func((pos, data)):
    if pos % 1000 == 0:
        print pos
    results['IDX'][pos] = pos
    results['CLUSZ'][pos] = data['CLUSZ'][0]
    results['LOSVD'][pos] = data['LOSVD'][0]
    results['LOSVDgmm'][pos] = data['LOSVDgmm'][0]
    results['R200'][pos] = data['R200'][0]
    results['MASS'][pos] = data['MASS'][0]

def find_indices(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i, _ in enumerate(smallArr):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds.append(sortedind[i1:i2])
        except IndexError:
            pass
        if i % 10000 ==0:
            print(i)

    return inds

def find_indices_multi(bigArr, smallArr, multi):
    from bisect import bisect_left, bisect_right
    from itertools import chain
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for j, sub in enumerate(multi):
        smallArr2 = smallArr[sub]
        tmp = []
        for i, _ in enumerate(smallArr2):
            i1 = bisect_left(sortedbigArr, smallArr2[i])
            i2 = bisect_right(sortedbigArr, smallArr2[i])
            try:
                tmp.append(sortedind[i1:i2])
            except IndexError:
                pass
        if j % 10000 == 0:
            print(j)

        inds.append(np.array(list(chain(*tmp))))
    return inds


def mkTruth():
    truthPath = './../data/buzzard_v1.0/allbands/truth/'
    # build truth database
    for i in range(20):
        with hdf.File(truthPath+'truth'+str(i).zfill(2)+'_Oii.hdf5', 'r') as f:
            dset = f[f.keys()[0]]
            print(dset.file) # print the loading file
            truth_part = dset['HALOID', 'RA', 'DEC', 'Z', 'Oii']
            try:
                truth = np.append(truth, truth_part)
            except NameError:
                truth = truth_part
    return truth

def mkHalo():
    haloPath = './../data/buzzard_v1.0/allbands/halos/'
    # build halo database
    for i in range(20):
        with hdf.File(haloPath+'halo'+str(i).zfill(2)+'.hdf5', 'r') as f:
            dset = f[f.keys()[0]]
            print(dset.file)
            halo_part = dset['id','upid', 'ra', 'dec', 'zspec', 'vrms', 'm200c']
            try:
                halo = np.append(halo, halo_part)
            except NameError:
                halo = halo_part

    return halo

if __name__ == "__main__":

    async_worker = AsyncFactory(worker, cb_func)
    halo = mkHalo()
    truth = mkTruth()

    mask = truth['g'] < 23
    truth = truth[mask]

    mask = (halo['m200c']/0.72 >= 1e13) & (halo['upid'] == -1)
    maskedHalo = halo[mask]
    hids, uniqueIdx = np.unique(maskedHalo['id'], return_index=True)

    # now we find all member halos of the most massive halos
    subHalos = find_indices(halo['upid'], hids)

    # now find all the corresponding galaxies
    gals = find_indices_multi(truth['HALOID'], halo['id'], subHalos)

    # make the results container
    results = np.zeros((len(subHalos),), dtype=[('IDX', '>i4'), ('HALOID',
        '>i8'), ('ZSPEC', '>f4'), ('VRMS', '>f4'), ('M200c', '>f4'), ('CLUSZ',
            '>f4'), ('LOSVD', '>f4'), ('LOSVDgmm', '>f4'), ('MASS', '>f4'),
            ('R200', '>f4'), ('NGAL', '>i4')])
    results['HALOID'] = hids

    print('do work')
    for i, SH in enumerate(subHalos):
        center = (maskedHalo['ra'][uniqueIdx[i]],
                maskedHalo['dec'][uniqueIdx[i]])
        if gals[i].size >= 5:
            async_worker.call(i, truth[gals[i]], center)
            # update results array
            results['NGAL'][i] = gals[i].size
            results['ZSPEC'][i] = maskedHalo['zspec'][uniqueIdx[i]]
            results['VRMS'][i] = maskedHalo['vrms'][uniqueIdx[i]]/np.sqrt(3)
            results['M200c'][i] = maskedHalo['m200c'][uniqueIdx[i]]/0.72
        else:
            results['IDX'][i] = i
            results['NGAL'][i] = gals[i].size
            results['ZSPEC'][i] = maskedHalo['zspec'][uniqueIdx[i]]
            results['VRMS'][i] = maskedHalo['vrms'][uniqueIdx[i]]/np.sqrt(3)
            results['M200c'][i] = maskedHalo['m200c'][uniqueIdx[i]]/0.72

    async_worker.wait()

    try:
        os.remove('result_targetedIdeal.hdf5')
    except OSError:
        pass
    with hdf.File('result_targetedIdeal.hdf5', 'w') as f:
        f['result_targetedIdeal'] = results
