import h5py as hdf
from functions import *
from sklearn.cluster import DBSCAN
from astLib import vec_astCalc
from astLib import astCoords
from numpy.lib import recfunctions as rfns
from calc_cluster_props import *
import pylab as pyl


# load RA/DEC/z data
f = hdf.File('../data/truth/Aardvark_v1.0c_truth_des_rotated.86.hdf5','r')
dset = f[list(f.keys())[0]]
truth = dset['RA', 'DEC', 'Z', 'HALOID']
f = hdf.File('../data/halos/Aardvark_v1.0_halos_r1_rotated.4.hdf5','r')
dset = f[list(f.keys())[0]]
halo = dset['HALOID', 'RA', 'DEC', 'Z', 'NGALS', 'M200']

# filter halos down
x = (halo['NGALS'] >= 5) & (halo['Z'] < 0.5) & (halo['M200'] >=1e13)
halo = halo[x]

# find the common halos
mask = pyl.in1d(truth['HALOID'], halo['HALOID'])
truth = truth[mask]

# find the haloids
haloids = pyl.intersect1d(halo['HALOID'], truth['HALOID'])

# find the indexes of the haloids in the halo list
inds = find_indices(halo['HALOID'], haloids)

# now we build the truth catalog with xx number of halos
# pick a few random halos
#randomHalos = pyl.random.choice(inds, 10)
randomHalos = pyl.array([42778,  1922,  2800,  5136,  9043, 42107, 42048, 37517,
    4399,  9184])

for i in randomHalos:
    x = pyl.where(truth['HALOID'] == halo['HALOID'][i])
    t_part = truth[x[0]]
    try:
        t = pyl.append(t, t_part)
    except NameError:
        t = t_part

    ra, dec, z = astCoords.eq2cart(halo['RA'][i], halo['DEC'][i],
            vec_astCalc.dm(halo['Z'][i]))
    print(ra, dec, z)
    #print halo['RA'][i], halo['DEC'][i], halo['Z'][i]
    print('members - ', len(x[0]))
    print('mass - ', halo['M200'][i])

print('---------')

# add some noise
noise = pyl.zeros(round(len(t)*0.6), dtype=t.dtype)
noise['RA'] = 340 + pyl.random(round(len(t)*0.6)) * 11
noise['DEC'] = -14 - pyl.random(round(len(t)*0.6)) * 8
noise['Z'] = pyl.random(round(len(t)*0.6)) * 0.5

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
rp = pyl.column_stack((X_,Y_,Z_))

# put it in a format it likes.
rp = pyl.array(rp.tolist())

# this does the finding
db = DBSCAN(eps=3, min_samples=6).fit(rp)
core_samples_mask = pyl.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
unique_labels = set(labels)

# how many clusters?
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

for k in unique_labels:
    if k == -1:
        break
    class_member_mask = (labels == k)
    xy = rp[class_member_mask & core_samples_mask]
    RA = pyl.mean(xy[:,0])
    DEC = pyl.mean(xy[:,1])
    Z = pyl.mean(xy[:,2])

    print('cluster - ', k, 'RA - ', RA, 'DEC - ', DEC, 'z - ', Z)
    print('members - ', len(xy))

    # update the data array with the recovered information
    t[class_member_mask & core_samples_mask] =\
    findClusterCenterRedshift(t[class_member_mask & core_samples_mask])

    t[class_member_mask & core_samples_mask] = findLOSV(t[class_member_mask &
        core_samples_mask])

    t['VD'][class_member_mask & core_samples_mask] =\
    calcVD_big(t['LOSV'][class_member_mask & core_samples_mask])

    t['MASS'][class_member_mask & core_samples_mask] =\
    calc_mass_Saro(t[class_member_mask & core_samples_mask])

### BELOW HERE IS ALL OF THE PLOTTING THAT WE DO ###




