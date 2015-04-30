import numpy as np
import h5py as hdf
from glob import glob
from calc_cluster_props import updateArray2
from halo_handler import find_indices, find_indices_bool

def update_result_file(hdfFile):
    # load the result file
    #f = hdf.File('out1204878.hdf5', 'r+')
    dset = hdfFile[hdfFile.keys()[0]]
    data = dset.value

    # update!
    data = updateArray2(data)

    # write out
    hdfFile['dset_appended'] = data

    hdfFile.flush()

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
        #f.close()
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

def fill_out_halo_info2(hdfFile):
    # load the result file
    #f = hdf.File('out1204878.hdf5', 'r+')
    dset = hdfFile[hdfFile.keys()[1]]
    data = dset.value

    # loads the halo files
    haloCat = mk_haloCatalog()
    # Get the unique haloids
    haloids = np.unique(data['HALOID'])
    # find the indexes for the halo information
    inds = find_indices_bool(haloCat['HALOID'], haloids)
    haloCat = haloCat[inds]

    # find the overlap. Will take a while.
    print 'find overlap'
    inds = find_indices(data['HALOID'], haloCat['HALOID'])

    print 'start loop'
    for idx, ind in enumerate(inds):
        data['CRA'][ind] = haloCat['RA'][idx]
        data['CDEC'][ind] = haloCat['DEC'][idx]
        data['CZ'][ind] = haloCat['Z'][idx]
        data['VRMS'][ind] = haloCat['VRMS'][idx]
        data['NGALS'][ind] = haloCat['NGALS'][idx]
        data['M200'][ind] = haloCat['M200'][idx]
        data['R200'][ind] = haloCat['R200'][idx]

        if idx % 10000 == 0:
            print idx

    hdfFile['dset_complete'] = data
    #f.close()
    hdfFile.flush()

    return data

if __name__ == '__main__':
    fill_out_halo_info2()
