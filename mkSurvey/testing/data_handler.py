import h5py as hdf
import numpy as np
#from numpy.lib import recfunctions
#data_dir = '../data/truth/'
data_dir = '/home/boada/scratch/hdf5/'

def load_tiles(tiles):
    data = []
    for t in tiles:
	t = t.replace('truth','')
        f = hdf.File(data_dir+'Aardvark_v1.0c_truth_des_rotated.'+ t +'.hdf5', 'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def mk_catalog(tiles):
    catalog = load_tiles(tiles)
    for dset in catalog:
        result_part = dset['ID', 'RA', 'DEC', 'Z', 'HALOID', 'RHALO', 'R200',
        'M200', 'NGALS']
        result_part = recfunctions.append_fields(result_part, 'OMAG',
                dset['OMAG'][:,1], usemask=False)
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part
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
        'M200', 'NGALS', 'OMAG']][selected]

    return result

