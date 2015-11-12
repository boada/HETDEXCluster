import numpy as np
import emcee
from astroML.datasets import generate_mu_z
from astroML.cosmology import Cosmology
import matplotlib.pyplot as plt
#from astLib import astCalc as aca


def compute_sigma_level(trace1, trace2, nbins=20):
    """From a set of traces, bin by number of standard deviations"""
    L, xbins, ybins = np.histogram2d(trace1, trace2, nbins)
    L[L == 0] = 1E-16
    logL = np.log(L)

    shape = L.shape
    L = L.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(L)[::-1]
    i_unsort = np.argsort(i_sort)

    L_cumsum = L[i_sort].cumsum()
    L_cumsum /= L_cumsum[-1]
    
    xbins = 0.5 * (xbins[1:] + xbins[:-1])
    ybins = 0.5 * (ybins[1:] + ybins[:-1])

    return xbins, ybins, L_cumsum[i_unsort].reshape(shape)


def plot_MCMC_trace(ax, xdata, ydata, trace, scatter=False, **kwargs):
    """Plot traces and contours"""
    xbins, ybins, sigma = compute_sigma_level(trace[0], trace[1])
    if scatter:
        ax.plot(trace[0], trace[1], 'ok', alpha=0.1, zorder=0)
    ax.contour(xbins, ybins, sigma.T, levels=[0.683, 0.955], **kwargs)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\beta$')
    
    
def plot_MCMC_model(ax, xdata, ydata, yerr, trace):
    """Plot the linear model and 2sigma contours"""
    #ax.plot(xdata, ydata, 'ok')
    ax.errorbar(xdata, ydata, yerr=yerr, fmt='ok')

    alpha, beta = trace[:2]
    xfit = np.linspace(xdata.min(), xdata.max(), 100)
    cosmo = Cosmology(h=np.median(trace[0]), omegaM=np.median(trace[1]), omegaL=np.median(trace[2]))
    yfit = np.asarray(map(cosmo.mu, xfit))
    yfit = yfit[:,np.newaxis]

    ax.plot(xfit, yfit, '-r')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def plot_MCMC_results(xdata, ydata, yerr, trace, colors='r'):
    """Plot both the trace and the model together"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    plot_MCMC_trace(ax[0], xdata, ydata, trace, True, colors=colors)
    plot_MCMC_model(ax[1], xdata, ydata, yerr, trace)


def calcMass(LOSVD, z, A1D=1082.9, alpha=0.3361):
    mass = 1e15 * (aca.H0 * aca.Ez(z)/100.) * (LOSVD/A1D)**(1/alpha)
    return mass

def log_prior(theta):
    H0, M0, L0 = theta
    if not 0.5 < H0 < 1:
        return -np.inf
    if not 0.05 < M0 < 0.75:
        return -np.inf
    if not 0.4 < L0 < 1.1:
        return -np.inf

    return 1

def log_likelihood(theta, z_sample, mu, dmu):
    H0, M0, L0 = theta
    print(theta)
    cosmo = Cosmology(omegaM=M0, omegaL=L0, h=H0)
    model = np.array(map(cosmo.mu, z_sample))
    return -0.5 * np.sum(np.log(2 * np.pi * dmu ** 2) + (mu - model) ** 2 /\
            dmu ** 2)

def log_posterior(theta, z_sample, mu, dmu):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, z_sample, mu, dmu)

# now we make some synethic data
z_sample, mu_sample, dmu = generate_mu_z(100, z0=0.3,
                                         dmu_0=0.05, dmu_1=0.004)

ndim = 3  # number of parameters in the model
nwalkers = 20  # number of MCMC walkers
nburn = 100  # "burn-in" period to let chains stabilize
nsteps = 500  # number of MCMC steps to take

# set theta near the maximum likelihood, with
np.random.seed()
starting_guesses = np.random.random((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[z_sample, mu_sample, dmu], threads=4)
sampler.run_mcmc(starting_guesses, nsteps)
print("done")

emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
plot_MCMC_results(z_sample, mu_sample, dmu, emcee_trace)


samples = sampler.chain[:, 2*nburn:, :].reshape((-1, ndim))

z_fit = np.linspace(0.04, 2, 100)
for h, m0, l0 in samples[np.random.randint(len(samples), size=50)]:
    cosmo = Cosmology(h=h, omegaM=m0, omegaL=l0)
    mu_fit = np.asarray(map(cosmo.mu, z_fit))
    plt.plot(z_fit, mu_fit, color="k", alpha=0.1)
cosmo = Cosmology(h=0.71, omegaM=0.27, omegaL=0.73)
mu_fit = np.asarray(map(cosmo.mu, z_fit))
plt.plot(z_fit, mu_fit , color="r", lw=2, alpha=0.8)
plt.errorbar(z_sample, mu_sample, yerr=dmu, fmt="ok")

m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
