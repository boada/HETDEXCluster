import numpy as np
import h5py as hdf
from glob import glob
from numpy.lib import recfunctions as rfns

def updateArray(data):
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    print 'update array...'
    newData = -np.ones(len(data))
    data = rfns.append_fields(data, ['CRA', 'CDEC', 'CZ', 'VRMS',
        'NGALS', 'M200', 'R200', 'CLUSZ', 'LOSV', 'LOSVD', 'MASS', 'Oii'],
        [newData, newData, newData, newData, newData, newData, newData,
            newData, newData, newData, newData, newData], dtypes='>f4',
        usemask=False)

    return data

def load_halos():
    ''' Loads all of the data sets. Doesn't actually load anything into
    memory just loads the top levels.

    '''

    print 'load halos...'
    data = []
    haloFiles = glob('../data/halos/choppedHaloFiles/*.hdf5')
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

    print 'make halo catalog...'
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
            print i

    return inds

def find_indices_single(bigArr, smallArr):
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
            inds.append(sortedind[i1:i2][0])
        except IndexError:
            pass
        if i % 10000 ==0:
            print i

    return inds


def find_indices_bool(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    It is important to remember that this returns a 1D list of True/False
    values which indicates which haloID is in the larger halo Catalog. This
    is different from the other find_indices which returns a multi-D list.

    '''

    inds = np.zeros(len(bigArr), dtype=bool)
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i, _ in enumerate(smallArr):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds[sortedind[i1:i2][0]] = True
        except IndexError:
            pass

    return inds

def fill_out_halo_info():
    # load the result file
    with hdf.File('out1204878.hdf5', 'r') as f:
        dset = f[f.keys()[0]]
        data = dset.value

        # update the data array to make room for the halo info
        data = updateArray(data)

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

    # now we write it out to a new file
    with hdf.File('out1204878_halo.hdf5', 'w') as f:
        f['dset_appended'] = data
        f.flush()

    return data

if __name__ == '__main__':
    fill_out_halo_info()


