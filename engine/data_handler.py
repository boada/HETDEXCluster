import h5py as hdf
import numpy as np
from numpy.lib import recfunctions as rfns

def find_tile(RAmin, DECmin, RAmax, DECmax, data=False):
    ''' Returns the name of the tile(s) that the current pointing is located
    inside of. The pointing is a box defined by the max/min of the RA/DEC.
    Throws an error if the pointing is wholely outside of the tiled region.

    '''

    if not len(data):
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

def mkTruth(i=-1, flatHMF=False):
    ''' Loads either all of the truth files and returns the data array or loads
    the individual file desired. To specify which file, call with an integer
    0-19 as an argument.

    '''
    truthPath = './../data/buzzard_v1.0/allbands/truth/'
    # build truth database
    if not i == -1:
        with hdf.File(truthPath+'truth'+str(i).zfill(2)+'_Oii.hdf5', 'r') as f:
            dset = f['truth'+str(i).zfill(2)+'_Oii']
            print(dset.file) # print the loading file
            if not flatHMF:
                truth_part = dset['HALOID', 'RA', 'DEC', 'Z', 'Oii']
            else:
                truth_part = dset['HALOID', 'RA', 'DEC', 'Z', 'Oii', 'VX',
                        'VY', 'VZ']
            truth_part = rfns.append_fields(truth_part, ['g','r'],
                    [dset['OMAG'][:,1], dset['OMAG'][:,2]], usemask=False)
            try:
                truth = np.append(truth, truth_part)
            except NameError:
                truth = truth_part
        return truth

    else:
        for i in range(20):
            with hdf.File(truthPath+'truth'+str(i).zfill(2)+'_Oii.hdf5', 'r')\
                as f:
                dset = f['truth'+str(i).zfill(2)+'_Oii']
                print(dset.file) # print the loading file
                if not flatHMF:
                    truth_part = dset['HALOID', 'RA', 'DEC', 'Z', 'Oii']
                else:
                    truth_part = dset['HALOID', 'RA', 'DEC', 'Z', 'Oii', 'VX',
                            'VY', 'VZ']
                truth_part = rfns.append_fields(truth_part, ['g','r'],
                        [dset['OMAG'][:,1], dset['OMAG'][:,2]],
                        usemask=False)
                try:
                    truth = np.append(truth, truth_part)
                except NameError:
                    truth = truth_part
        return truth

def mkHalo():
    ''' Loads *ALL* of the halo information and returns the data array. '''

    haloPath = './../data/buzzard_v1.0/allbands/halos/'
    # build halo database
    for i in range(20):
        with hdf.File(haloPath+'halo'+str(i).zfill(2)+'.hdf5', 'r') as f:
            dset = f[f.keys()[0]]
            print(dset.file)
            halo_part = dset['id','upid', 'ra', 'dec', 'zspec', 'vrms',
                    'm200c', 'rvir']
            try:
                halo = np.append(halo, halo_part)
            except NameError:
                halo = halo_part

    return halo

def mkQs(i=-1):
    ''' Loads either all of the truth files and returns the data array or loads
    the individual file desired. To specify which file, call with an integer
    0-19 as an argument.

    '''
    truthPath = './../data/buzzard_v1.0/allbands/truth/'
    # build truth database
    if not i == -1:
        with hdf.File(truthPath+'truth'+str(i).zfill(2)+'_Oii.hdf5', 'r') as f:
            dset = f['Q']
            print(dset.file) # print the loading file
            q_part = dset.value
            try:
                q = np.append(q, q_part)
            except NameError:
                q = q_part
        return q

    else:
        for i in range(20):
            with hdf.File(truthPath+'truth'+str(i).zfill(2)+'_Oii.hdf5', 'r')\
                as f:
                dset = f['Q']
                print(dset.file) # print the loading file
                q = dset.value
                try:
                    q = np.append(q, q_part)
                except NameError:
                    q = q_part
        return q

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
    #result = catalog[['HID', 'RA', 'DEC', 'Z', 'HALOID', 'g', 'r']][selected]
    result = catalog[selected]

    return result

