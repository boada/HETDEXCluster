import numpy as np
import h5py as hdf
from functions import *
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from astLib import vec_astCalc
from astLib import astCoords
from astLib import astStats
from numpy.lib import recfunctions as rfns
from calc_cluster_props import *

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
    sigmaLocal = astStats.biweightScale(localData['LOSV'], tuningConstant=9.0)
    vLocal = astStats.biweightLocation(localData['LOSV'], tuningConstant=6.0)

    Nnn = int(np.sqrt(len(localData)))

    delta2 = (Nnn + 1)/ sigma**2 * ((vLocal - v)**2 + (sigmaLocal - sigma)**2)

    return delta2

def DStest():
    """ Computes the delta deviation for an entire cluster. See Dressler et al.
    1988 for a complete description. A result close to the number of galaxies
    present in the cluster indicates no substructure.

    """

    dsm = DSmetric()


# load RA/DEC/z data
f = hdf.File('../data/truth/Aardvark_v1.0c_truth_des_rotated.86.hdf5','r')
dset = f[f.keys()[0]]
truth = dset['RA', 'DEC', 'Z', 'HALOID']
f = hdf.File('../data/halos/Aardvark_v1.0_halos_r1_rotated.4.hdf5','r')
dset = f[f.keys()[0]]
halo = dset['HALOID', 'RA', 'DEC', 'Z', 'NGALS', 'M200']

# filter halos down
x = (halo['NGALS'] >= 5) & (halo['Z'] < 0.5) & (halo['M200'] >=1e13)
halo = halo[x]

# find the common halos
mask = np.in1d(truth['HALOID'], halo['HALOID'])
truth = truth[mask]

# find the haloids
haloids = np.intersect1d(halo['HALOID'], truth['HALOID'])

# find the indexes of the haloids in the halo list
inds = find_indices(halo['HALOID'], haloids)

# now we build the truth catalog with xx number of halos
# pick a few random halos
#randomHalos = np.random.choice(inds, 10)
randomHalos = np.array([42778,  1922,  2800,  5136,  9043, 42107, 42048, 37517,
    4399,  9184])

for i in randomHalos:
    x = np.where(truth['HALOID'] == halo['HALOID'][i])
    t_part = truth[x[0]]
    try:
        t = np.append(t, t_part)
    except NameError:
        t = t_part

    ra, dec, z = astCoords.eq2cart(halo['RA'][i], halo['DEC'][i],
            vec_astCalc.dm(halo['Z'][i]))
    print ra, dec, z
    #print halo['RA'][i], halo['DEC'][i], halo['Z'][i]
    print 'members - ', len(x[0])
    print 'mass - ', halo['M200'][i]

print '---------'

# add some noise
noise = np.zeros(round(len(t)*0.6), dtype=t.dtype)
noise['RA'] = 340 + np.random.rand(round(len(t)*0.6)) * 11
noise['DEC'] = -14 - np.random.rand(round(len(t)*0.6)) * 8
noise['Z'] = np.random.rand(round(len(t)*0.6)) * 0.5

# add the noise to the data array
t = rfns.stack_arrays((t,noise), usemask=False)

# add extra fields to the array
t = updateArray2(t)

# add the mass information
for i in randomHalos:
    x = t['HALOID'] == halo['HALOID'][i]
    t['M200'][x] = halo['M200'][i]

# convert to physical units
X_,Y_,Z_ = astCoords.eq2cart(t['RA'], t['DEC'], vec_astCalc.dm(t['Z']))

# now we stack all the parts to pass to the cluster finder
rp = np.column_stack((X_,Y_,Z_))

# put it in a format it likes.
rp = np.array(rp.tolist())

# this does the finding
db = DBSCAN(eps=3, min_samples=6).fit(rp)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
unique_labels = set(labels)

# how many clusters?
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

# this updates the t matrix
for k in unique_labels:
    if k == -1:
        break
    class_member_mask = (labels == k)
    xy = rp[class_member_mask & core_samples_mask]
    RA = np.mean(xy[:,0])
    DEC = np.mean(xy[:,1])
    Z = np.mean(xy[:,2])

    print 'cluster - ', k, 'RA - ', RA, 'DEC - ', DEC, 'z - ', Z
    print 'members - ', len(xy)

    # update the data array with the recovered information
    t[class_member_mask & core_samples_mask] =\
    findClusterCenterRedshift(t[class_member_mask & core_samples_mask])

    t[class_member_mask & core_samples_mask] = findLOSV(t[class_member_mask &
        core_samples_mask])

    t['VD'][class_member_mask & core_samples_mask] =\
    calcVD_big(t['LOSV'][class_member_mask & core_samples_mask])

    t['MASS'][class_member_mask & core_samples_mask] =\
    calc_mass_Saro(t[class_member_mask & core_samples_mask])

# now we do all of the DS test stuff

for k in unique_labels:
    if k == -1:
        break

    # find the cluster we are going to work with
    class_member_mask = (labels == k)
    c = rp[class_member_mask & core_samples_mask]
    # find the nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=int(np.sqrt(len(c))),
        algorithm='ball_tree').fit(c)
    distances, indices = nbrs.kneighbors(c)

    t2 = t[class_member_mask & core_samples_mask]

    ds2 = [(t2['RA'][inds[0]], t2['DEC'][inds[0]], DSmetric(t2[inds],
        t2['LOSV'].mean(), t2['VD'][0])) for inds in
            indices]

    ds = np.sqrt(ds2)

    delta = np.sum(ds)
    print delta/len(c)




