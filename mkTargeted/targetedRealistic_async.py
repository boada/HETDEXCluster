from multiprocessing import Pool
import h5py as hdf
import numpy as np
from calc_cluster_props import *
from data_handler import mkTruth, mkHalo
import os
import sys

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

def worker(pos, data, center, tZ):
    #print "IN:PID: %d \t Value: %d" % (os.getpid(), pos)
    data = updateArray(data)
    #data = findClusterRedshift(data)
    data['CLUSZ'] = tZ
    data = findSeperationSpatial(data, center)
    data = findLOSV(data)
    # make initial cuts
    mask = abs(data['LOSV']) < 5000
    data = data[mask]
    while True:
        try:
            if size == data.size:
                break
        except NameError:
            pass

        size = data.size
        #print 'size', data.size

        #data = rejectInterlopers(data)
        try:
            x = shifty_gapper(data['SEP'], data['Z'], tZ, ngap=15, glimit=500)
            data = data[x]
        except:
            break
        #data = findLOSVD(data)
        data = findLOSVDgmm(data)
        data['LOSVD'] = data['LOSVDgmm']

        data = findR200(data)
        mask = data['SEP'] < data['R200'][0]
        data = data[mask]

        data = findClusterRedshift(data)
        data = findSeperationSpatial(data, center)

    #data = findLOSVDgmm(data)
    data = calc_mass_Evrard(data, A1D = 1177, alpha = 0.364)
    #print "OUT:PID: %d \t Value: %d" % (os.getpid(), pos)
    return pos, data

def cb_func((pos, data)):
    if pos % 1000 == 0:
        print pos
    results['IDX'][pos] = pos
    results['CLUSZ'][pos] = data['CLUSZ'][0]
    results['LOSVD'][pos] = data['LOSVD'][0]
    results['LOSVDgmm'][pos] = data['LOSVDgmm'][0]
    results['MASS'][pos] = data['MASS'][0]
    results['R200'][pos] = data['R200'][0]
    results['NGAL'][pos] = data.size

if __name__ == "__main__":
    async_worker = AsyncFactory(worker, cb_func)
    halo = mkHalo()
    truth = mkTruth()

    mask = truth['g'] < 23.
    truth = truth[mask]

    mask = (halo['m200c']/0.72 >= 1e13) & (halo['upid'] == -1)
    maskedHalo = halo[mask]
    hids, uniqueIdx = np.unique(maskedHalo['id'], return_index=True)

    # limit the cases
    print sys.argv[1], sys.argv[2]
    if int(sys.argv[1]) == 0:
        hids = hids[:int(sys.argv[2])]
        uniqueIdx = uniqueIdx[:int(sys.argv[2])]
    elif int(sys.argv[2]) == 9:
        hids = hids[int(sys.argv[1]):]
        uniqueIdx = uniqueIdx[int(sys.argv[1]):]
    else:
        hids = hids[int(sys.argv[1]):int(sys.argv[2])]
        uniqueIdx = uniqueIdx[int(sys.argv[1]):int(sys.argv[2])]

    # make the results container
    results = np.zeros((hids.size,), dtype=[('IDX', '>i4'), ('HALOID',
        '>i8'), ('ZSPEC', '>f4'), ('VRMS', '>f4'), ('M200c', '>f4'), ('RVIR',
        '>f4'), ('CLUSZ', '>f4'), ('LOSVD', '>f4'), ('LOSVDgmm', '>f4'),
        ('MASS', '>f4'), ('R200', '>f4'), ('NGAL', '>i4')])
    results['HALOID'] = hids
    # now we have to make some initial cuts and then make final spatial cuts
    for i, SH in enumerate(hids):
        center = (maskedHalo['ra'][uniqueIdx[i]],
            maskedHalo['dec'][uniqueIdx[i]])
        raMask = (center[0] - 0.5 < truth['RA']) & (truth['RA'] < center[0] + 0.5)
        decMask = (center[1] - 0.5 < truth['DEC']) & (truth['DEC'] < center[1] +
                0.5)
        async_worker.call(i, truth[raMask & decMask], center,
                maskedHalo['zspec'][uniqueIdx[i]])

        results['ZSPEC'][i] = maskedHalo['zspec'][uniqueIdx[i]]
        results['VRMS'][i] = maskedHalo['vrms'][uniqueIdx[i]]/np.sqrt(3)
        results['M200c'][i] = maskedHalo['m200c'][uniqueIdx[i]]/0.72
        results['RVIR'][i] = maskedHalo['rvir'][uniqueIdx[i]]/0.72

    async_worker.wait()

    with hdf.File('result_targetedRealistic'+str(os.environ['LSB_JOBID'])+'.hdf5',
            'w') as f:
        f['result_targetedRealistic'] = results
