import h5py as hdf
import numpy as np
#data_dir = '../data/truth/'
data_dir = '/home/boada/scratch/truth/'

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
            f = hdf.File(data_dir+'Aardvark_v1.0_halos_r1_rotated.'+t+'.hdf5',
                    'r')
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
            f = hdf.File(data_dir+'Aardvark_v1.0c_truth_des_rotated.'+ t\
                    +'.hdf5', 'r')
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
    result = catalog[['ID', 'RA', 'DEC', 'Z', 'HALOID', 'RHALO', 'R200',
        'M200', 'OMAG']][selected]

    return result

