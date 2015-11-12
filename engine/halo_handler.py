import h5py as hdf
import numpy as np

data_dir = '/home/boada/scratch/halos/'
#data_dir= '/Users/steven/Projects/desCluster/data/halos/'

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

def find_indices_bool(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

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

def findRADECmaxmin(ra, dec):
    ''' Figures out the max and min RA/DEC for the observed catalog. Use this
    info to figure out which halo tiles we need to look into to get the rest of
    the infomation. Returns RA/DEC max RA/DEC min.

    '''

    if ra.max() > 300. and ra.min() < 100:
        x = np.where(ra < 200)
        y = np.where(ra > 200)
        return ra[x].max(), dec.max(), ra[y].min(), dec.min()
    else:
        return ra.max(), dec.max(), ra.min(), dec.min()

def findHalotile(RAmin, DECmin, RAmax, DECmax, data=False):
    ''' Returns the name of the tile(s) that the current pointing is located
    inside of. The pointing is a box defined by the max/min of the RA/DEC.
    Throws an error if the pointing is wholely outside of the tiled region.

    '''

    if not len(data):
        data = fix_haloTiles()
        #data = np.genfromtxt('tiles.txt', names=True, dtype=None)
    else:
        pass
    # Find all of the tiles that don't overlap with the box defined. We'll take
    # the inverse of that below to find all of the tiles we want.
    # left/right
    leftRight = (RAmin > data['RAmax'] ) | (RAmax < data['RAmin'])
    # top/bottom
    topBottom = (DECmin > data['DECmax']) | (DECmax < data['DECmin'])

    tile = np.intersect1d(data['name'][~leftRight],data['name'][~topBottom])

    if len(tile):
        return tile
    else:
        raise ValueError('Out of RA/DEC bounds!')

def fix_haloTiles():
    ''' fixes the tile catalog for the tiles that overlap zero. Now the tiles
    will go from RAmin -> 360 and 0 -> RAmax. Still should only return one tile
    if only one tile is covering the data point.

    '''

    data = np.genfromtxt('halos.txt', names=True, dtype=None)
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
