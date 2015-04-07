import h5py as hdf
import numpy as np

#data_dir = '/home/boada/scratch/'
data_dir= '/Users/steven/Projects/desCluster/data/halos/'

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

def mk_haloCatalog(tiles):
    ''' This does all the actual loading of the data. Only selects out the
    columns that we think we'll need to work with.

    '''

    catalog = load_halos(tiles)
    for dset in catalog:
        result_part = dset['HALOID', 'RA', 'DEC', 'Z', 'VRMS', 'NGALS', 'M200',
                'R200']
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    return result

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

def findHaloTile(RAmax, DECmax, RAmin, DECmin, data=False):
    ''' Figure out what halo tiles we are going to need to load. Works with the
    max and min RA/DEC positions to figure that out.

    '''

    if not len(data):
        data = fix_haloTiles()
    else:
        pass
    tileDEC = (data['DECmax'] > DECmin) & (DECmin > data['DECmin'])
    tileRA = (data['RAmax'] > RAmin) & (RAmin > data['RAmin'])
    tile = np.intersect1d(data['name'][tileRA], data['name'][tileDEC])

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
