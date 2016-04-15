from __future__ import division
import h5py as hdf
from scipy import stats
import numpy as np
import emcee

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

def log_prior(theta):
    m, b, s = theta
    if s >= 0:
        return 0
    return -np.inf

def log_likelihood(theta, x, y, yerr):
    m, b, s = theta
    model = m*x + b
    sigma = s**2 + yerr**2

    return -0.5 * np.sum(np.log(2 *np.pi*sigma) + (y-model)**2 / sigma)

def log_probfn(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, x, y, yerr)

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
    target = dset['HALOID', 'M200c']

with hdf.File('./targetedRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err']

# mask out the values with failed ML masses
mask = (target_mlMasses['ML_pred_3d'] != 0)
target = target[mask]
target_mlMasses = target_mlMasses[mask]

with hdf.File('./surveyCompleteRealistic.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey = dset['HALOID', 'M200c']

with hdf.File('./surveyCompleteRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err']

# mask out the values with failed ML masses
mask = (survey_mlMasses['ML_pred_3d'] != 0)
survey = survey[mask]
survey_mlMasses = survey_mlMasses[mask]

scatter = 0.05

for d, m, in zip([target, survey], [target_mlMasses, survey_mlMasses]):
    bins = np.arange(10,150,20)

    # add the noise to the true masses-- 0.25 dex at the moment
    m_obs = stats.norm(np.log10(d['M200c']), scatter).rvs(d.size)
    # use the noisy masses to calculate an observed lambda
    lam_obs = mklambda(truth)(m_obs)

    # setup the data
    x_obs = np.log(lam_obs/40.)
    y_obs = np.log(10**m['ML_pred_3d']/1e14)

    yerr = np.zeros_like(y_obs)

    # Set up the sampler.
    nwalkers, ndim = 100, 3
    p0 = np.random.random((nwalkers, ndim))
    #p0 = np.append(truth, x_obs)
    #p0 = [p0 + np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probfn,
                                    args=[x_obs, y_obs, yerr])

    # Burn in.
    print("Burning in.")
    p0, lnprob0, state = sampler.run_mcmc(p0, 500)
    sampler.reset()

    # Sample.
    print("Sampling.")
    sampler.run_mcmc(p0, 1000)

    # Print results.
    samples = sampler.flatchain
    print("m = {0} ± {1}".format(np.mean(samples[:, 0]), np.std(samples[:, 0])))
    print("b = {0} ± {1}".format(np.mean(samples[:, 1]), np.std(samples[:, 1])))





