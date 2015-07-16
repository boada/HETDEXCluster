import h5py as hdf
import numpy as np
from numpy.lib import recfunctions
#data_dir = '../data/truth/'
data_dir = '/home/boada/scratch/truth/'

def find_tile(RAmin, DECmin, RAmax, DECmax, data=False):
    ''' Returns the name of the tile(s) that the current pointing is located
    inside of. The pointing is a box defined by the max/min of the RA/DEC.
    Throws an error if the pointing is wholely outside of the tiled region.

    '''

    if data is False:
        data = fix_tiles()
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

def fix_tiles():
    ''' fixes the tile catalog for the tiles that overlap zero. Now the tiles
    will go from RAmin -> 360 and 0 -> RAmax. Still should only return one tile
    if only one tile is covering the data point.

    '''

    data = np.genfromtxt('buzzard_truth.txt', names=True, dtype=None)
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

def load_tiles(tiles):
    data = []
    for t in tiles:
	t = t.replace('truth','')
        f = hdf.File(data_dir+'Buzzard-highres_galaxies_shmatch.'+ t +'.hdf5',
                'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def mk_catalog(tiles):
    catalog = load_tiles(tiles)
    for dset in catalog:
        print dset.file # The file it is loading
        result_part = dset['ID', 'RA', 'DEC', 'Z', 'HALOID']
        result_part = recfunctions.append_fields(result_part, ['g','r'],
                [dset['OMAG'][:,0], dset['OMAG'][:,1]], usemask=False)
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    print 'done loading'
    return result

def apply_mask(ramin, decmin, ramax, decmax, catalog):
    if ramin < ramax:
        x = (ramax > catalog['RA']) & (catalog['RA'] > ramin)
    else:
        x1 = (ramax > catalog['RA']) & (catalog['RA'] > 0)
        x2 = (360 > catalog['RA']) & (catalog['RA'] > ramin)
        # Merge the results together
        x = x1 | x2
    y = (decmax > catalog['DEC']) & (catalog['DEC'] > decmin)

    selected = x & y
    result = catalog[['ID', 'RA', 'DEC', 'Z', 'HALOID', 'g', 'r']][selected]

    return result

