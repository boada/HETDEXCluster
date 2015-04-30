import h5py as hdf
import numpy as np

def fix_tiles(tiles=None):
    ''' fixes the tile catalog for the tiles that overlap zero. Now the tiles
    will go from RAmin -> 360 and 0 -> RAmax. Still should only return one tile
    if only one tile is covering the data point.

    '''

    if tiles == None:
        data = np.genfromtxt('tiles.txt', names=True, dtype=None)
    else:
        data = np.genfromtxt(tiles, names=True, dtype=None)

    # Find the tiles that overlap zero.
    x = np.where(data['RAmax'] < data['RAmin'])
    d2 = data[x]
    d2['RAmax'] = 360.
    d3 = data[x]
    d3['RAmin'] = 0.0
    data = np.delete(data, x)
    data = np.append(data, d2)
    data = np.append(data, d3)

    return data

def mk_haloCatalog(tiles):
    ''' This does all the actual loading of the data. Only selects out the
    columns that we think we'll need to work with.

    '''

    data_dir = '/home/boada/scratch/halos/'

    def load_halos(tiles):
        ''' Loads all of the data sets. Doesn't actually load anything into
        memory just loads the top levels.

        '''

        data = []
        for t in tiles:
            t = t.replace('halos','')
            f = hdf.File(data_dir+'Aardvark_v1.0_halos_r1_rotated.'+t+'.hdf5', 'r')
            dset = f[f.keys()[0]]
            data.append(dset)
        return data

    catalog = load_halos(tiles)
    for dset in catalog:
        result_part = dset['HALOID', 'VRMS', 'Z', 'M200', 'R200', 'NGALS']
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    return result

def mk_catalog(tiles):

    data_dir = '/home/boada/scratch/truth/'

    def load_tiles(tiles):
        data = []
        for t in tiles:
            t = t.replace('truth','')
            f = hdf.File(data_dir+'Aardvark_v1.0c_truth_des_rotated.'+ t +'.hdf5', 'r')
            dset = f[f.keys()[0]]
            data.append(dset)
        return data

    catalog = load_tiles(tiles)
    for dset in catalog:
        print dset.file # The file it is loading
        result_part = dset['HALOID', 'Z', 'RHALO', 'R200']
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    print 'done loading'
    return result

def find_indices(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i in range(len(smallArr)):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds.append(sortedind[i1:i2][0])
        except IndexError:
            pass

    return inds
