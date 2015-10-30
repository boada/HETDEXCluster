import numpy as np
import emcee
from astLib import astCalc as aca

def calcMass(LOSVD, z, A1D=1082.9, alpha=0.3361):
    mass = 1e15 * (aca.H0 * aca.Ez(z)/100.) * (LOSVD/A1D)**(1/alpha)
    return mass

def log_prior(theta):
    H0, M0, L0, sigma = theta
    if sigma < 0:
        return -np.inf  # log(0)
    if H0 < 60 or M0 < 0 or L0 < 0:
        return -np.inf
    if not M0 + L0 ==1:
        return -np.inf
    return 1

def log_likelihood(theta, y, z, LOSVD):
    H0, M0, L0, sigma = theta
    aca.H0 = H0
    aca.OMEGA_M0 = M0
    aca.OMEGA_L0 = L0
    y_model = np.log10(np.array(map(calcMass, LOSVD, z)))
    return -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (y - y_model) ** 2 /\
            sigma ** 2)

def log_posterior(theta, y, z, LOSVD):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, y, z, LOSVD)

# choose true parameters
aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

# now we make some synethic data
N = 50
LOSVD = 200 + 1000*np.random.rand(N)
z = 0.1 + 0.4 * np.random.rand(N)
mTrue = np.array(map(calcMass, LOSVD, z))
mObs = mTrue + mTrue * np.random.normal(0,0.5,size=N)
m_err = 0.1 + 0.5 * np.random.rand(N)
mObs += m_err * np.random.randn(N)

ndim = 4  # number of parameters in the model
nwalkers = 50  # number of MCMC walkers
nburn = 1000  # "burn-in" period to let chains stabilize
nsteps = 2000  # number of MCMC steps to take

# set theta near the maximum likelihood, with
np.random.seed(0)
starting_guesses = np.random.random((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[mTrue,
    z, LOSVD])
sampler.run_mcmc(starting_guesses, nsteps)
print("done")

emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
#plot_MCMC_results(xdata, ydata, emcee_trace)
