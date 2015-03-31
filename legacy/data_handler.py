import h5py as hdf
import numpy as np
from numpy.lib import recfunctions
data_dir = './data/truth/'
#data_dir = '/home/boada/scratch/hdf5/'

def load_tiles(tiles):
    data = []
    for t in tiles:
	t = t.replace('truth','')
        f = hdf.File(data_dir+'Aardvark_v1.0c_truth_des_rotated.'+ t +'.hdf5', 'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def apply_mask(ramin, decmin, ramax, decmax, tiles):
    catalog = load_tiles(tiles)
    print tiles
    for tileID, name in enumerate(tiles):
        print tileID, name
        dset = catalog[tileID]
        if ramin < ramax:
            x = (ramax > dset['RA']) & (dset['RA'] > ramin)
        else:
            x1 = (ramax > dset['RA']) & (dset['RA'] > 0)
            x2 = (360 > dset['RA']) & (dset['RA'] > ramin)
            # Merge the results together
            x = x1 | x2
        y = (decmax > dset['DEC']) & (dset['DEC'] > decmin)

        selected = x & y

        if len(tiles) > 1:
            result_part = dset['ID', 'RA', 'DEC', 'Z', 'HALOID', 'RHALO',
                    'R200', 'M200', 'NGALS'][selected]
            result_part = recfunctions.append_fields(result_part, 'OMAG',
                    dset['OMAG'][:,1][selected], usemask=False)
            try:
                result = np.append(result, result_part)
            except NameError:
		result = result_part
        else:
            result = dset['ID', 'RA', 'DEC', 'Z', 'HALOID', 'RHALO', 'R200',
                'M200', 'NGALS'][selected]
            result = recfunctions.append_fields(result, 'OMAG',
                    dset['OMAG'][:,1][selected], usemask=False)

    return result

