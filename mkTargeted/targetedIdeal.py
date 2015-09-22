import multiprocessing
import numpy as np
import h5py as hdf
from calc_cluster_props import updateArray, findClusterRedshift,\
    findSeperationSpatial, findLOSV
from os import environ

def mp_worker((gals, center)):
    # get the galaxy properties
    data = truth[gals]
    # make room for other calculations
    data = updateArray(data)
    data = findClusterRedshift(data)
    data = findSeperationSpatial(data, center)
    data = findLOSV(data)


def mp_worker_wrapper(args):
    return mp_worker(*args)

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

    #mask = halo['m200c']/0.72 >= 1e13
    #return halo[mask]
    return halo

if __name__ == "__main__":
    p = multiprocessing.Pool(multiprocessing.cpu_count())

    halo = mkHalo()
    truth = mkTruth()

    mask = (halo['m200c'] >= 1e13) & (halo['upid'] == -1) & (halo['zspec'] <\
            0.2)
    maskedHalo = halo[mask]
    # find all of the halos which are already the most massive
    #mask = maskedHalo['upid'] == -1
    hids = maskedHalo['id']
    # now we find all member halos of the most massive halos
    subHalos = find_indices(halo['upid'], hids[:10])

    # this is where the loop happens
    gals = np.ravel(find_indices(truth['HALOID'], halo['id'][subHalos[0]]))
    props = truth[gals]



