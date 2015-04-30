import multiprocessing
import numpy as np
import h5py as hdf
from glob import glob
from astLib.astCoords import shiftRADec
# my functions
from halo_handler import *
import calc_cluster_props

def mp_worker(result, clusRA, clusDEC, clusZ):
    ''' This is going to do all of the determinations for the individual
    cluster properties. Thinking we'll feed it a copy of the full
    catalog/result and have it return the properties of the single cluster.

    '''

    upper = shiftRADec(clusRA, clusDEC, 300, 300)
    lower = shiftRADec(clusRA, clusDEC, -300, -300)

    x1 = (upper[0] > result['RA']) & (result['RA'] > lower[0])
    x2 = (upper[1] > result['DEC']) & (result['DEC'] > lower[1])
    selected = x1 & x2

    r = result[selected]

    # make a velocity space cut.
    y = (abs(calc_cluster_props.findLOSV(r, clusZ)) < 10000)
    r = r[y]

    if len(r)  <6:
        return -1, -1, -1
    else:
        r = calc_cluster_props.updateArray(r)
        r = calc_cluster_props.findClusterCenterRedshift(r)
        r = calc_cluster_props.findLOSV(r)
        if  len(r) >=15:
            r = calc_cluster_props.findSeperationSpatial(r, (clusRA, clusDEC),
                    unit='Mpc')
            r = calc_cluster_props.rejectInterlopers(r)
        else:
            print 'group!'
            r = calc_cluster_props.findSeperationSpatial(r, (clusRA, clusDEC),
                    unit='arcsecond')
            r = calc_cluster_props.rejectInterlopers_group(r)

        # fewer than 6 members left... return -1
        if len(r) < 6:
            return -1, -1 , -1
        else:
            vd, m200, r200 = calc_cluster_props.calc_mass(r)
            #return calc_cluster_props.findClusterCenterRedshift(r)
            return vd, m200, r200

def mp_worker_wrapper(args):
    ''' wrapper for the worker function. '''

    return mp_worker(*args)

def mp_handler(catalog, result):
    from itertools import repeat, izip, imap

    # find the individual haloids that we'll be working with
    haloids = np.unique(result['HALOID'])
    # finds the indices of the halos that we found in the observations in the
    # halo catalog file.
    inds = find_indices(catalog['HALOID'], haloids)

    result = map(mp_worker, repeat(result, len(inds)), catalog['RA'][inds],
        catalog['DEC'][inds], catalog['Z'][inds])

#    p = multiprocessing.Pool()
#    result = p.map(mp_worker_wrapper, izip(repeat(result, len(inds)),
#        catalog['RA'][inds], catalog['DEC'][inds], catalog['Z'][inds]))
#    p.close()
#    p.join()

    return result


if __name__ == "__main__":
    # load result file
    f = hdf.File('test.hdf5', 'r')
    #f = hdf.File('out0.176032713212.hdf5', 'r')
    dset = f[f.keys()[0]]
    result = dset.value

    # have to filter by redshift and r-mag
    x = (0.47 >result['Z']) & (22. > result['OMAG'])
    result = result[x]

    # figure out the max/min RA/DEC to know which tiles to load
    coords = findRADECmaxmin(result['RA'], result['DEC'])
    data = fix_haloTiles()
    tiles = findHaloTile(coords[0], coords[1], coords[2], coords[3], data=data)

    print tiles
    # load the halo tiles into a catalog that we can work with
    catalog = mk_haloCatalog(tiles)
    # filter the catalog to only have big clusters and the right redshift
    x = (catalog['M200'] >= 1e13) & (0.47 > catalog['Z'])
    catalog = catalog[x]

    done = mp_handler(catalog, result)
