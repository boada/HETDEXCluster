import numpy as np
import h5py as hdf
from numpy.lib import recfunctions
from bisect import bisect_left, bisect_right
import multiprocessing
from itertools import repeat
data_dir = '/home/boada/scratch/truth/'

def fix_truthTiles():
    ''' fixes the tile catalog for the tiles that overlap zero. Now the tiles
    will go from RAmin -> 360 and 0 -> RAmax. Still should only return one tile
    if only one tile is covering the data point.

    '''

    data = np.genfromtxt('../mkSurvey/truth.txt', names=True, dtype=None)
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

def findTruthTile(RAmax, DECmax, RAmin, DECmin):
    ''' Figure out what halo tiles we are going to need to load. Works with the
    max and min RA/DEC positions to figure that out.

    '''

    data = fix_truthTiles()
    # left/right
    leftRight = (RAmin > data['RAmax'] ) | (RAmax < data['RAmin'])
    # top/bottom
    topBottom = (DECmin > data['DECmax']) | (DECmax < data['DECmin'])

    tile = np.intersect1d(data['name'][~leftRight],data['name'][~topBottom])

    if len(tile):
        return tile
    else:
        raise ValueError('Out of RA/DEC bounds!')

def load_tiles(tiles):
    data = []
    for t in tiles:
        t = t.replace('truth','')
        f = hdf.File(data_dir+'Aardvark_v1.0c_truth_des_rotated.'+ t +'.hdf5',
                'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def mk_catalog(tiles):
    catalog = load_tiles(tiles)
    for dset in catalog:
        print dset.file # The file it is loading
        result_part = dset['RA', 'DEC', 'Z', 'HALOID']
        result_part = recfunctions.append_fields(result_part, ['g','r'],
                [dset['OMAG'][:,0], dset['OMAG'][:,1]], usemask=False)
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
    print 'done loading'
    return result

def find_indices(bigArr, smallArr):
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
            if len(sortedind[i1:i2]) !=0:
                inds.append(sortedind[i1:i2])
        except IndexError:
            pass

    return inds

def worker((hids, tiles)):
    print tiles
    truth = mk_catalog([tiles])
    # try to make smaller
    mask = truth['Z'] < 0.5
    truth = truth[mask]

    # now we find what we are looking for, could take a while
    print 'finding halos', multiprocessing.current_process().name
    halos = find_indices(truth['HALOID'], hids)

    print len(halos)
    if not len(halos):
        print 'NO HALOS FOUND!'
        return []

    for h in halos:
        result_part = truth[h]
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part

    return result

def worker_wrapper(*args):
    return worker(*args)

def start_process():
    print 'Starting', multiprocessing.current_process().name

# setup the data we selected with the HETDEX mask
with hdf.File('out1204878_complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    selectedData = dset['HALOID', 'M200', 'Z', 'CRA', 'CDEC']

mask = (selectedData['M200']/0.72 >= 1e13) & (selectedData['Z'] < 0.5)
selectedData = selectedData[mask]
hids = np.unique(selectedData['HALOID'])

# find the truth tiles we think we need
tiles = findTruthTile(selectedData['CRA'].max(), selectedData['CDEC'].max(),
        selectedData['CRA'].min(), selectedData['CDEC'].min())
print tiles
# and load that data

p = multiprocessing.Pool(10, maxtasksperchild=2, initializer=start_process)

result = p.map(worker_wrapper, zip(repeat(hids), tiles))

p.close()
p.join()

for r in result:
    if len(r) > 0:
        try:
            final = np.append(final, r)
        except NameError:
            final = r
    else:
        'Halos?'
with hdf.File('out1204878_allGalaxies.hdf5', 'w') as f:
    f['ds_complete'] = final
    f.flush()


