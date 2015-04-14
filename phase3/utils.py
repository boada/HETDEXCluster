import numpy as np
import h5py as hdf
from glob import glob
from calc_cluster_props import updateArray2
from halo_handler import find_indices_bool

def update_result_file():
    # load the result file
    f = hdf.File('out1204878.hdf5', 'r+')
    dset = f[f.keys()[0]]
    data = dset.value

    # update!
    data = updateArray2(data)

    # write out
    f['dset_appended'] = data

    f.close()

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

def fill_out_halo_info():
    # load the result file
    f = hdf.File('out1204878.hdf5', 'r+')
    dset = f[f.keys()[1]]
    data = dset.value

    # Get the unique haloids
    haloids = np.unique(data['HALOID'])

    # loads the halo files
    haloCat = mk_haloCatalog()

    # find the indexes for the halo information
    inds = find_indices_bool(haloCat['HALOID'], haloids)

    haloCat = haloCat[inds]

    haloids_inds = np.argsort(data['HALOID'])
    haloCat_inds = np.argsort(haloCat['HALOID'])

    b = []
    i = 0
    for hci in haloCat_inds:
        look = haloCat['HALOID'][hci]
        while look != data['HALOID'][haloids_inds[i]]:
            if look < data['HALOID'][haloids_inds[i]]:
                print 'overshot!'
                raise ValueError
            i+=1
        while 1:
            try:
                if look == data['HALOID'][haloids_inds[i]]:
                    data['CRA'][haloids_inds[i]] = haloCat['HALOID'][hci]
                    data['CDEC'][haloids_inds[i]] = haloCat['DEC'][hci]
                    data['CZ'][haloids_inds[i]] = haloCat['Z'][hci]
                    data['VRMS'][haloids_inds[i]] = haloCat['VRMS'][hci]
                    data['NGALS'][haloids_inds[i]] = haloCat['NGALS'][hci]
                    data['M200'][haloids_inds[i]] = haloCat['M200'][hci]
                    data['R200'][haloids_inds[i]] = haloCat['R200'][hci]
                    i+=1
                else:
                    #print i, look, data['HALOID'][haloids_inds[i]], 'between'
                    b.append(i)
                    break
            except IndexError:
                break
        #if i % 5000 == 0:

#    f['dset_complete'] = data
#    f.close()
    print len(b)

    return data

if __name__ == '__main__':
    fill_out_halo_info()
