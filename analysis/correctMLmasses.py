import h5py as hdf
import numpy as np
import sys


def main(f1, f2, f3):
    ### biased ###
    ###############
    with hdf.File(f1, 'r') as f:
        dset = f[f.keys()[0]]
        data = dset.value
    # filter bad values
    mask = (data['ML_pred_1d'] != 0)
    data = data[mask]

    ### corrections ###
    ###################
    with hdf.File(f2, 'r') as f:
        dset = f[f.keys()[0]]
        data_bias = dset.value

    # these are the mass bins we used for the corrections
    bins = np.arange(11.5, 16, 0.1)

    # make the results container
    results = np.copy(data)

    results['M200c'] = np.log10(data['M200c'])

    indx = np.digitize(np.log10(data['MASS']), bins)
    for i in range(1, bins.size):
        results['MASS'][indx==i] =\
        np.log10(data['MASS'][indx==i]) - data_bias['powerlaw_bias'][i-1]

    indx = np.digitize(data['ML_pred_1d'], bins)
    for i in range(1, bins.size):
        results['ML_pred_1d'][indx==i] =\
        data['ML_pred_1d'][indx==i] - data_bias['ML_bias_1d'][i-1]

    indx = np.digitize(data['ML_pred_3d'], bins)
    for i in range(1, bins.size):
        results['ML_pred_2d'][indx==i] =\
        data['ML_pred_2d'][indx==i] - data_bias['ML_bias_2d'][i-1]

    indx = np.digitize(data['ML_pred_3d'], bins)
    for i in range(1, bins.size):
        results['ML_pred_3d'][indx==i] =\
        data['ML_pred_3d'][indx==i] - data_bias['ML_bias_3d'][i-1]

    with hdf.File(f3, 'w') as f:
        f['corrected'] = results


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
