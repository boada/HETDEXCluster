import h5py as hdf
import numpy as np
from sklearn import mixture
from addHaloInfo import find_indices
from astLib import astStats

#with hdf.File('out1204878_hetdex.hdf5', 'r') as f:
with hdf.File('out1204878_complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

mask = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.5)
data = data[mask]

hids = np.unique(data['HALOID'])
halos = np.array(find_indices(data['HALOID'], hids))

newLOSVDS = np.zeros([6000, 5])

for i,h in enumerate(halos[:6000]):
    if len(h) < 10:
        pass
    else:

        LOSV = data['LOSV'][h]

        lowest_bic = np.infty
        bic = []
        n_components_range = range(1, 4)
        cv_types = ['spherical', 'tied', 'diag', 'full']
        for cv_type in cv_types:
            for n_components in n_components_range:
                # Fit a mixture of Gaussians with EM
                gmm = mixture.GMM(n_components=n_components,
                        covariance_type=cv_type)
                gmm.fit(LOSV)
                bic.append(gmm.bic(LOSV))
                if bic[-1] < lowest_bic:
                    lowest_bic = bic[-1]
                    best_gmm = gmm

        # figure things out
        covars = best_gmm.covars_.ravel()
        weights = best_gmm.weights_.ravel()
        means = best_gmm.means_.ravel()
        wmeans = np.sum(weights*means)

        parts = weights * ((means - wmeans)**2 + covars)
        newLOSVD = np.sqrt(np.sum(parts))


#        parts = best_gmm.weights_.ravel()**2 * best_gmm.covars_.ravel()
#        newLOSVD = np.sum(np.sqrt(parts))

        resample = best_gmm.sample(1000)
        newnewLOSVD = astStats.biweightScale_test(resample, tuningConstant=9.)
        # now we resample and then see
        #dx = np.linspace(LOSV.min()-100,LOSV.max()+100,1000)
        #logprob, responsibilities = best_gmm.eval(dx)
        #pdf = np.exp(logprob)

#        print best_gmm.weights_.ravel()
#        print data['VRMS'][h][0]/np.sqrt(3), data['LOSVD'][h][0], newLOSVD, newnewLOSVD

        newLOSVDS[i] = [h.shape[0],data['VRMS'][h][0]/np.sqrt(3),
                data['LOSVD'][h][0], newLOSVD, newnewLOSVD]

