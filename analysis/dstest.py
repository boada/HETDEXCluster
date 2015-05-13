import numpy as np
from sklearn.neighbors import NearestNeighbors
from calc_cluster_props import *
from astLib import astStats

def DSmetric(localData, v, sigma):
    """ Calculates the delta squared values as given in equation 1 of Dressler
    et al. 1988. This uses both the position and velocities to give a measure
    of substructure by identifying galaxies with do not follow the cluster
    velocity distribution.

    @type localData: list
    @param localData: The velocities of the target galaxy and its N nearest
        neighbors.
    @type v: float
    @param v: The global mean velocity of the cluster
    @type sigma: float
    @param sigma: The global velocity dispersion of the cluster
    @rtype: list
    @return: a list of delta squared values

    """

    if 3 <= localData.shape[0] < 15:
        sigmaLocal = astStats.gapperEstimator(localData)
    elif 15 <= localData.shape[0]:
        sigmaLocal = astStats.biweightScale(localData, tuningConstant=9.0)

    try:
        vLocal = astStats.biweightLocation(localData, tuningConstant=6.0)
        #vLocal = localData.mean()
        if vLocal == None:
            vLocal = localData.mean()
    except ZeroDivisionError:
        print 'exception'
        vLocal = localData.mean()

    Nnn = localData.shape[0]

    delta2 = (Nnn + 1)/ sigma**2 * ((vLocal - v)**2 + (sigmaLocal - sigma)**2)

    return delta2

def DStest(data, LOSV, LOSVD, method='shuffle', shuffles=1e4):
    """ Computes the delta deviation for an entire cluster. See Dressler et al.
    1988 for a complete description. A result close to the number of galaxies
    present in the cluster indicates no substructure.

    @type data: array
    @param data: The 2D (3D) array of RA, DEC, (z) galaxy positions.
    @type LOSV: array
    @param LOSV: Array of line of sight velocities (LOSV) for each galaxy in
        cluster
    @type LOSVD: float
    @param LOSVD: Global line of sight velocity dispersion for the cluster.
    @type method: string
    @param method: The method used for the determination of substructure.
        Should either be 'shuffle' (default) or 'threshold'. 'Shuffle' shuffles
        the velocities and computes the probability of substructure.
        'threshold' uses a simple threshold and is generally considered not as
        reliable as the shuffle method.
    @type shuffles: float
    @param shuffles: The number of times to shuffle the velocities.
    @rtype: float
    @rparam: The significance of substructure. For method='shuffle' this is the
        probability of substructure, and for method='threshold' a value greater
        than 1 indicates the presence of substructure.

    """

    try:
        data.shape[1]
    except IndexError:
        raise IndexError('data must be at least 2D')

    nbrs = NearestNeighbors(n_neighbors=int(np.sqrt(len(data))),
        algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)

    ds2 = [DSmetric(LOSV[inds], LOSV.mean(), LOSVD) for inds in indices]
    ds = np.sqrt(ds2)
    delta = np.sum(ds)

    if method == 'threshold':
        return delta/data.shape[0]

    elif method == 'shuffle':
        deltaShuffled = np.zeros(shuffles)
        for i in range(int(shuffles)):
            np.random.shuffle(LOSV)
            ds2 = [DSmetric(LOSV[inds], LOSV.mean(), LOSVD) for inds in
                indices]
            ds = np.sqrt(ds2)
            deltaShuffled[i] = np.sum(ds)

            if i+1%1000 ==0:
                print i

        mask = deltaShuffled > delta

        return deltaShuffled[mask].shape[0]/float(shuffles)

    else:
        raise NameError("method must be either 'threshold' or 'shuffle'")

def findLOSV(data, CLUSZ):
    ''' Finds the line of sight velocity for each of the galaxies and puts it
    in the LOSV column of the data array.

    '''
    c = 2.99e5 # speed of light in km/s
    losv = c * (data - CLUSZ)/(1 + CLUSZ)
    return losv

def main():
    N = 100
    N2 = 1


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

    data = np.column_stack([ra,dec])

    s = DStest(data, LOSV, LOSVD, shuffles=3000)
    return s
if __name__ == "__main__":
    main()
