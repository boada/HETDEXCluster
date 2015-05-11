import numpy as np
from sklearn.neighbors import NearestNeighbors
from calc_cluster_props import *
from astLib import astStats
from addHaloInfo import find_indices

def DSmetric(localData, v, sigma):
    """ Calculates the delta squared values as given in equation 1 of Dressler
    et al. 1988. This uses both the position and velocities to give a measure
    of substructure by identifying galaxies with do not follow the cluster
    velocity distribution.

    @type x: list
    @param x: x (RA) position, accepts a list of RA values in decimal degrees
    @type y: list
    @param y: y (DEC) position, accepts a list of DEC values in decimal degrees
    @type v: list
    @param v: velocities associated with the galaxies
    @rtype: list
    @return: a list of delta squared values

    """
    if 3 <= len(localData) < 15:
        sigmaLocal = astStats.gapperEstimator(localData)
    elif 15 <= len(h):
        sigmaLocal = astStats.biweightScale(localData, tuningConstant=9.0)

    try:
        vLocal = astStats.biweightLocation(localData, tuningConstant=6.0)
    except ZeroDivisionError:
        vLocal = localData.mean()

    Nnn = int(np.sqrt(len(localData)))+1

    delta2 = (Nnn + 1)/ sigma**2 * ((vLocal - v)**2 + (sigmaLocal - sigma)**2)

    return delta2

#with hdf.File('out1204878_hetdex.hdf5', 'r') as f:
with hdf.File('out1204878_complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

mask = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.5)
data = data[mask]

hids = np.unique(data['HALOID'])
halos = np.array(find_indices(data['HALOID'], hids))

blah = []
for h in halos:
    if len(h) < 10:
        pass
    else:
        cluster = np.column_stack([data['RA'][h], data['DEC'][h], data['LOSV'][h]])
        nbrs = NearestNeighbors(n_neighbors=int(np.sqrt(len(cluster)))+1,
            algorithm='ball_tree').fit(cluster)
        distances, indices = nbrs.kneighbors(cluster)

        ds2 = np.array([DSmetric(cluster[:,2][inds], cluster[:,2].mean(),
            data['LOSVD'][h][0]) for inds in indices])

        ds = np.sqrt(ds2)

        delta = np.sum(ds)
        blah.append([len(h),delta/len(cluster)])
#    print delta/len(cluster)
