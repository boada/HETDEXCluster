from halo_handler import *
from calc_cluster_props import *
import h5py as hdf
import numpy as np
import multiprocessing
from os import environ

def make_work(truth, halo, inds):
    for i in range(len(inds)):
    #for i in range(10):
        x = np.where(truth['HALOID'] == halo['HALOID'][inds[i]])
        t = truth[x[0]]
        a = np.where(t['RHALO'] < t['R200'])
        t = t[a[0]]
        t = rfns.append_fields(t, ['CVRMS', 'CR200', 'CM200'],
            [halo['VRMS'][inds[i]]/np.sqrt(3), halo['R200'][inds[i]],
            halo['M200'][inds[i]]], usemask=False)
        yield t

def mp_handler(truth, halo, inds):

    p = multiprocessing.Pool()
    result = p.map(mp_worker, make_work(truth, halo, inds))
    p.close()
    p.join()

    #result = map(mp_worker, make_work(truth, halo, inds))

    return result

def mp_worker(t):

    lt = len(t)
    if lt >= 5:
        id = t['HALOID'][0]
        #print id
        t = updateArray(t)
        t = findClusterCenterRedshift(t)
        t = findLOSV(t)

        vd_uc = np.std(t['LOSV'])
        vd_c = np.std(t['LOSV'], ddof=1)
        vd_g = np.std(t['LOSV'], ddof=1.5)
        vd_bi = calcVD_big(t['LOSV'])
        vd_gap = calcVD_small(t['LOSV'])
        vd_t = calcVD_test(t['LOSV'])

        return (id, lt,vd_uc, vd_c, vd_g, vd_bi, vd_gap, vd_t, t['CVRMS'][0],
                t['CM200'][0], t['CLUSZ'][0])

    else:
        return -1, -1, -1, -1, -1, -1, -1, -1, -1 , -1 , -1

if __name__ == "__main__":
    base = '/home/boada/scratch'

    f = hdf.File(base+'/truth/Aardvark_v1.0c_truth_des_rotated.61.hdf5')
    dset = f[f.keys()[0]]
    truth = dset['HALOID', 'Z', 'RHALO', 'R200', 'M200']
    f.close()

    f = hdf.File(base+'/halos/Aardvark_v1.0_halos_r1_rotated.4.hdf5')
    dset = f[f.keys()[0]]
    halo = dset['HALOID', 'Z', 'VRMS', 'R200', 'M200']
    f.close()

    # Only big clusters:
    x = (halo['M200'] >= 1e13)
    halo = halo[x]

    haloids = np.unique(truth['HALOID'])
    inds = find_indices(halo['HALOID'], haloids)

    print 'do work'
    result = mp_handler(truth, halo, inds)

    # write it all to a file
    print 'ID number is', environ['LSB_JOBID']
    f = hdf.File('stat_testing'+str(environ['LSB_JOBID'])+'.hdf5', 'w')
    f['dset'] = result
    f.close()
    print 'done'

