import numpy as np
from sklearn.neighbors import NearestNeighbors
from calc_cluster_props import *
from astLib import astStats

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
    sigmaLocal = astStats.biweightScale(localData, tuningConstant=9.0)
    vLocal = astStats.biweightLocation(localData, tuningConstant=6.0)

    Nnn = int(np.sqrt(len(localData)))

    delta2 = (Nnn + 1)/ sigma**2 * ((vLocal - v)**2 + (sigmaLocal - sigma)**2)

    return delta2

def DStest():
    """ Computes the delta deviation for an entire cluster. See Dressler et al.
    1988 for a complete description. A result close to the number of galaxies
    present in the cluster indicates no substructure.

    """

    dsm = DSmetric()

def findLOSV(data, CLUSZ):
    ''' Finds the line of sight velocity for each of the galaxies and puts it
    in the LOSV column of the data array.

    '''
    c = 2.99e5 # speed of light in km/s
    losv = c * (data - CLUSZ)/(1 + CLUSZ)
    return losv


N = 100
N2 = 20

# make some data
ra1 = 0 + np.random.rand(N) * 0.5
dec1 = -1 - np.random.rand(N) * 0.5
z1 = 0.2 + np.random.rand(N) * 0.01

# add another little cluster
ra2 = 0.4 + np.random.rand(N2) * 0.1
dec2 = -1.1 - np.random.rand(N2) * 0.1
z2 = 0.21 + np.random.rand(N2) * 0.01

ra = np.append(ra1, ra2)
dec = np.append(dec1, dec2)
z = np.append(z1, z2)



#calculate things we would have normally
CLUSZ = astStats.biweightLocation(z, tuningConstant=6.0)
LOSV = findLOSV(z, CLUSZ)
LOSVD = astStats.biweightScale_test(LOSV, tuningConstant=9.0)


data = np.column_stack([ra,dec,z])


nbrs = NearestNeighbors(n_neighbors=int(np.sqrt(len(data))),
    algorithm='ball_tree').fit(data)
distances, indices = nbrs.kneighbors(data)

ds2 = np.array([(ra[inds[0]], dec[inds[0]], DSmetric(LOSV[inds],
    LOSV.mean(), LOSVD)) for inds in indices])

ds = np.sqrt(ds2[:,-1])

delta = np.sum(ds)
print delta/(N+N2)
