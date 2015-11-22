import numpy as np
import emcee
from astLib import astCalc as aca
from astLib import astStats as ast

aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def findLOSVD(LOSV):
    if LOSV.size >=15:
        return ast.biweightScale_test(LOSV, tuningConstant=9.0)
    else:
        return ast.gapperEstimator(LOSV)

def calcMass(vd, A1D = 1082.9, alpha=0.3361):
    avgz = 0.0
    return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

def log_prior(theta, LOSV):
    sigma, mu = theta
    if not 0 < sigma < 1200:
        return -np.inf
    if not LOSV.min() < mu < LOSV.max():
        return -np.inf

    return 1

def log_likelihood(theta, LOSV, LOSV_err):
    sigma, mu = theta
    #print(theta)

    # break long equation into three parts
    a = -0.5 * np.sum(np.log(LOSV_err**2 + sigma**2))
    b = -0.5 * np.sum((LOSV - mu)**2/(LOSV_err**2 + sigma**2))
    c = -1. * LOSV.size/2. * np.log(2*np.pi)

    return a +b +c

def log_posterior(theta, LOSV, LOSV_err):
    lp = log_prior(theta, LOSV)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, LOSV, LOSV_err)

def mkLOSVS(vd, ngals):
    return np.random.normal(loc=0, scale=vd, size=int(ngals))

def plot_training(mass, sigma, bins=41):
    gridx = np.linspace(mass.min(), mass.max(), bins+1)
    gridy = np.linspace(sigma.min(), sigma.max(), bins+1)
    H, xbins, ybins = np.histogram2d(mass, sigma, bins=[gridx, gridy])

    return H, xbins, ybins

N_clusters = 5000

# make the truth
LOSVDS = np.random.uniform(100, 1200, N_clusters)
NGALS = np.random.randint(5, 250, N_clusters)
MASS = np.array(map(calcMass, LOSVDS))

# Now we recover the LOSVDs and this becomes our training sample
sigmas_obs = np.zeros_like(LOSVDS)
for i, (vd, ngals) in enumerate(zip(LOSVDS, NGALS)):
    LOSV = mkLOSVS(vd, ngals)
    # right now they have zero error
    LOSV_err = np.zeros_like(LOSV)
    sigmas_obs[i] = findLOSVD(LOSV)

# Now we try to recover, and predict
expected_mass = np.zeros(N_clusters)
for j, (vd, ngals) in enumerate(zip(LOSVDS[:10], NGALS[:10])):
    LOSV = mkLOSVS(vd, ngals)
    # right now they have zero error
    LOSV_err = np.zeros_like(LOSV)

    # here is the MCMC stuff
    ndim = 2  # number of parameters in the model
    nwalkers = 40  # number of MCMC walkers
    nburn = 50  # "burn-in" period to let chains stabilize
    nsteps = 300  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    m = np.random.normal(np.mean(LOSV), scale=1, size=(nwalkers))
    s = np.random.normal(np.std(LOSV), scale=1, size=(nwalkers))
    starting_guesses = np.vstack([s,m]).T

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[LOSV,
        LOSV_err], threads=4)
    sampler.run_mcmc(starting_guesses, nsteps)

    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
#    sigma_rec, mean_rec = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
#                             zip(*np.percentile(samples, [16, 50, 84],
#                                                axis=0)))

    # using the distribution of sigmas we make P(s) and ds
    Ps, ds = np.histogram(np.log10(samples[:,0]), bins=50, density=True)
    # and get new new sigmas from the histogram2d
    centers = (ds[1:] + ds[:-1])/2.0

    # now we have to make P(mt|s) -- from the training sample
    # will need to do this a bunch of times for each sigma in centers
    H, xbins, ybins = plot_training(np.log10(MASS), np.log10(sigmas_obs))
    H = H.T
    iy  = np.digitize(centers, ybins) #change this line
    # now we make P(mt) using P(mt|s) and P(s)ds from above
    Pmt = np.zeros(xbins.size -1)
    Psds = Ps*np.diff(ds)
    for i in range(iy.size):
        try:
            Pmt_s = H[iy[i]] / H[iy[i]].sum()
        except IndexError:
            Pmt_s = np.zeros(xbins.size -1)
        Pmt += Pmt_s * Psds[i]

    # now we can calculate the expected mass
    centers = (xbins[:-1] + xbins[1:])/2.
    norm = np.sum(Pmt * np.diff(xbins))
    M_expect = np.sum(centers * Pmt * np.diff(xbins))/norm

    print np.log10(MASS[j]),M_expect

