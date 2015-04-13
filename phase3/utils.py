import h5py as hdf
from calc_cluster_props import updateArray2

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
