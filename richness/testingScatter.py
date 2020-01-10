
from scipy import stats
import numpy as np


def mklogMass(theta):
    ''' Makes Log10 Mass from a given richness following the simet2016
    mass-richness relationship.

    '''
    m0, alpha, lambda0 = theta
    return lambda x: m0 + alpha * np.log10(x / lambda0)


def mklambda(theta):
    ''' Makes Lambda given a Log10 Mass. If mass is not Log10 it won't work
    correctly.

    '''

    m0, alpha, lambda0 = theta
    return lambda y: lambda0 * (10**y / 10**m0)**(1 / alpha)


# normalization, power-law index, and lambda0 of the lambda-mass relation
truth = 14.344, 1.33, 40  # simet2016
# truth = 14.191, 1.31, 30 # farahi2016

recovery = []
for i in range(1000):
    ### Fake Data for Testing ###
    N = 1000
    x_err, y_err = 0.5, 0.1
    #m_true = np.random.uniform(13,15.5, N)
    m_true = np.linspace(13, 15.5, N)
    m_pert = stats.norm(m_true, x_err).rvs(N)
    m_pred = stats.norm(m_true, y_err).rvs(N)  # unbiased data
    m_pred_bias = stats.norm(m_true, y_err).rvs(N) + np.linspace(0.75, 0,
                                                                N)  # biased data
    lam = mklambda(truth)(m_true)

    lam_bins = np.arange(10, 140, 10)
    idx = np.digitize(lam, lam_bins)
    for i in range(1, lam_bins.size):
        recovery.append(np.std(m_pred[idx == i]))
        #print(np.std(m_pred[idx == i]))
