import numpy as np
import h5py as hdf
from glob import glob
from calc_cluster_props import updateArray2
from halo_handler import find_indices_bool
import multiprocessing

def update_result_file():
    # load the result file
    f = hdf.File('out1204878.hdf5', 'r+')
    dset = f[f.keys()[0]]
    data = dset.value

    # update!
    data = updateArray2(data)

    # write out
    f['dset_appended'] = data

    f.close()

def load_halos():
    ''' Loads all of the data sets. Doesn't actually load anything into
    memory just loads the top levels.

    '''

    data = []
    haloFiles = glob('./haloFiles/*.hdf5')
    for f in haloFiles:
        f = hdf.File(f, 'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def mk_haloCatalog():
    ''' This does all the actual loading of the data. Only selects out the
    columns that we think we'll need to work with.

    '''

    catalog = load_halos()
    for dset in catalog:
        print dset.file
        result_part = dset['HALOID', 'RA', 'DEC', 'Z', 'VRMS', 'NGALS', 'M200',
                'R200']
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    return result

def mp_worker(idx, hid):
    ii = np.where(haloids == hid)[0]
    data['CRA'][ii] = haloCat['RA'][idx]
    data['CDEC'][ii] = haloCat['DEC'][idx]
    data['CZ'][ii] = haloCat['Z'][idx]
    data['VRMS'][ii] = haloCat['VRMS'][idx]
    data['NGALS'][ii] = haloCat['NGALS'][idx]
    data['M200'][ii] = haloCat['M200'][idx]
    data['R200'][ii] = haloCat['R200'][idx]
    if idx % 1000 == 0:
        print idx

    return 0

def mp_worker_wrapper(args):
    return mp_worker(*args)


def fillOutInfo_multi():
    # load the result file
    f = hdf.File('out1204878.hdf5', 'r+')
    dset = f[f.keys()[1]]
    myglobals.data = dset.value

    # Get the unique haloids
    haloids = np.unique(myglobals.data['HALOID'])

    # loads the halo files
    myglobals.haloCat = mk_haloCatalog()

    # find the indexes for the halo information
    inds = find_indices_bool(myglobals.haloCat['HALOID'], haloids)

    myglobals.haloCat = myglobals.haloCat[inds]

    #haloids = data['HALOID']

    p = multiprocessing.Pool(2)
    result = p.map(mp_worker_wrapper, enumerate(haloCata['HALOID']))
    p.join()
    p.close()

    f['dset_complete'] = myglobals.data
    f.close()

def fill_out_halo_info():
    # load the result file
    f = hdf.File('out1204878.hdf5', 'r+')
    dset = f[f.keys()[1]]
    data = dset.value

    # Get the unique haloids
    haloids = np.unique(data['HALOID'])

    # loads the halo files
    haloCat = mk_haloCatalog()

    # find the indexes for the halo information
    inds = find_indices_bool(haloCat['HALOID'], haloids)

    haloCat = haloCat[inds]

    haloids = data['HALOID']

    for idx, hid in enumerate(haloCat['HALOID']):
        ii = np.where(haloids == hid)[0]
        data['CRA'][ii] = haloCat['RA'][idx]
        data['CDEC'][ii] = haloCat['DEC'][idx]
        data['CZ'][ii] = haloCat['Z'][idx]
        data['VRMS'][ii] = haloCat['VRMS'][idx]
        data['NGALS'][ii] = haloCat['NGALS'][idx]
        data['M200'][ii] = haloCat['M200'][idx]
        data['R200'][ii] = haloCat['R200'][idx]
        if idx % 1000 == 0:
            print idx
        haloids = np.delete(haloids, ii)

    f['dset_complete'] = data
    f.close()

if __name__ == '__main__':
    #fill_out_halo_info()
    fillOutInfo_multi()
