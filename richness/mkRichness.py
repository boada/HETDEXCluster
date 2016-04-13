from __future__ import division
import h5py as hdf
from scipy import stats
import numpy as np


def mklogMass(theta):
    ''' Makes Log10 Mass from a given richness following the simet2016
    mass-richness relationship.

    '''
    m0, alpha, lambda0 = theta
    return lambda x: m0 + alpha*np.log10(x/lambda0)

def mklambda(theta):
    ''' Makes Lambda given a Log10 Mass. If mass is not Log10 it won't work
    correctly.

    '''

    m0, alpha, lambda0 = theta
    return lambda y: lambda0*(10**y/10**m0)**(1/alpha)

# normalization, power-law index, and lambda0 of the lambda-mass relation
truth = 14.344, 1.33, 40

### Fake Data for Testing ###
#N = 100
#x_true = np.linspace(1,250, N)
#y_true = mklogMass(truth)(x_true)

#x_err, y_err = 0, 0.250
#x_obs = stats.norm(x_true, x_err).rvs(N)
#y_obs = stats.norm(y_true, y_err).rvs(N)

with hdf.File('./result_targetedRealistic.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset['HALOID', 'M200c']

m_obs = stats.norm(np.log10(data['M200c']), 0.25).rvs(data.size)
lam_obs = mklambda(truth)(m_obs)


with hdf.File('./targetedRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    mlMasses = dset['HALOID', 'ML_pred_3d']




