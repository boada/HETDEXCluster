
import h5py as hdf
from scipy import stats
import numpy as np
import emcee


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


def log_prior(theta):
    m, b, s = theta
    if s >= 0:
        return 0
    return -np.inf


def log_likelihood(theta, x, y, yerr):
    m, b, s = theta
    model = m * x + b
    sigma = s**2 + yerr**2

    return -0.5 * np.sum(np.log(2 * np.pi * sigma) + (y - model)**2 / sigma)


def log_probfn(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, x, y, yerr)

# normalization, power-law index, and lambda0 of the lambda-mass relation
truth = 14.344, 1.33, 40  # simet2016
# truth = 14.191, 1.31, 30 # farahi2016
# truth = 14.2042, 1.1655, 30 # me

with hdf.File('./targetedRealistic_MLmasses_corrected.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    target_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err', 'M200c']

# mask out the values with failed ML masses
mask = (target_mlMasses['ML_pred_3d'] != 0)
target_mlMasses = target_mlMasses[mask]

with hdf.File('./surveyCompleteRealistic_MLmasses_corrected.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    survey_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err', 'M200c']

# mask out the values with failed ML masses
mask = (survey_mlMasses['ML_pred_3d'] != 0)
survey_mlMasses = survey_mlMasses[mask]

scatter = 0.25

for m in [target_mlMasses, survey_mlMasses]:
    # add the noise to the true masses-- 0.25 dex at the moment
    m_obs = stats.norm(m['M200c'], scatter).rvs(m.size)
    # use the noisy masses to calculate an observed lambda
    lam_obs = mklambda(truth)(m_obs)
    mask = (10 <= lam_obs) & (lam_obs < 130)

    # setup the data
    x_obs = np.log10(lam_obs)[mask]
    y_obs = m['ML_pred_3d'][mask]
    print(np.std(m['M200c'][mask] - y_obs))
    print(np.std(m_obs[mask] - y_obs))

    yerr = m['ML_pred_3d_err'][mask]
    #yerr = np.zeros_like(y_obs)
    #yerr = np.ones_like(y_obs) * scatter

    # Set up the sampler.
    nwalkers, ndim = 100, 3
    p0 = np.random.random((nwalkers, ndim))
    #p0 = np.append(truth, x_obs)
    #p0 = [p0 + np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers,
                                    ndim,
                                    log_probfn,
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
    print('m = {0} +/- {1}'.format(
        np.median(samples[:, 0]), np.std(samples[:, 0])))
    print('b = {0} +/- {1}'.format(
        np.median(samples[:, 1]), np.std(samples[:, 1])))
    print('s = {0} +/- {1}'.format(
        np.median(samples[:, 2]), np.std(samples[:, 2])))
