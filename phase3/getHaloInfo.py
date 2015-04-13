import h5py as hdf
import numpy as np
from halo_handler import find_indices_bool
from glob import glob
from calc_cluster_props import updateArray2

def load_halos():
    ''' Loads all of the data sets. Doesn't actually load anything into
    memory just loads the top levels.

    '''

    data = []
    haloFiles = glob('./haloFiles/*.hdf5')
    for f in haloFiles:
        f = hdf.File(f, 'r')
        dset = f[f.keys()[0]]
        data.append(dset)
    return data

def mk_haloCatalog():
    ''' This does all the actual loading of the data. Only selects out the
    columns that we think we'll need to work with.

    '''

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



# load the result file
f = hdf.File('out1204878.hdf5', 'r+')
dset = f[f.keys()[0]]




# Get the unique haloids
haloids = np.unique(dset['HALOID'])

# loads the halo files
haloCat = mk_haloCatalog()


inds = find_indices_bool(haloCat['HALOID'], haloids)

data = dset.value




