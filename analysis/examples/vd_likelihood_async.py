import numpy as np
import emcee
from multiprocessing import Pool

class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool(maxtasksperchild=10)

    def call(self,*args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()


def worker(LOSV, LOSV_err):
    # here is the MCMC stuff
    ndim = 2  # number of parameters in the model
    nwalkers = 40  # number of MCMC walkers
    nburn = 100  #  "burn-in" period to let chains stabilize
    nsteps = 500  # number of MCMC steps to take

    np.random.seed()
    #starting_guesses = np.random.random((nwalkers, ndim))
    m = np.random.normal(np.mean(LOSV), scale=1, size=(nwalkers))
    s = np.random.normal(np.std(LOSV), scale=1, size=(nwalkers))
    starting_guesses = np.vstack([s,m]).T

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[LOSV,
        LOSV_err])
    sampler.run_mcmc(starting_guesses, nsteps)

    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    sigma_rec, mean_rec = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0)))

    return mean_rec, sigma_rec

def cb_func((mean_rec, sigma_rec)):
    resultsMean.append(mean_rec)
    resultsSigma.append(sigma_rec)

def log_prior(theta, LOSV):
    sigma, mu = theta
    if  sigma < 0:
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

if __name__ == "__main__":
    # multiprocessing stuff
    async_worker = AsyncFactory(worker, cb_func)

    mean, sigma = 0, 400

    Ntrials = 10
    resultsMean = []
    resultsSigma = []
    for i in range(Ntrials):
        # make fake data
        LOSV = np.random.normal(loc=mean, scale=sigma, size=100)
        # start with zero errors
        LOSV_err = np.zeros_like(LOSV)
        #LOSV_err = 0.1+0.5*np.random.rand(LOSV.size)

        async_worker.call(LOSV, LOSV_err)

    async_worker.wait()

    resultsMean = np.array(resultsMean)
    resultsSigma = np.array(resultsSigma)

    meanDown = resultsMean[:,0] - resultsMean[:,2]
    meanUp = resultsMean[:,0] + resultsMean[:,1]
    sigmaDown = resultsSigma[:,0] - resultsSigma[:,2]
    sigmaUp = resultsSigma[:,0] + resultsSigma[:,1]

    meanMask = (meanDown < mean) & (mean < meanUp)
    sigmaMask = (sigmaDown < sigma) & (sigma < sigmaUp)

    print('Mean percentage', resultsMean[meanMask][:,0].size/float(Ntrials))
    print('Sigma percentage', resultsMean[sigmaMask][:,0].size/float(Ntrials))

