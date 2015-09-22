import h5py as hdf
import pylab as pyl
from scipy.integrate import quad
from astLib import astCalc

data_dir = '../data/buzzard_v1.0/allbands/halos/'

# load all of the data
for i in range(20):
    print i
    with hdf.File(data_dir+'halo'+str(i).zfill(2)+'.hdf5', 'r') as f:
        dset = f[f.keys()[0]]
        result_part = dset['m200c', 'zspec']
        try:
            result = pyl.append(result, result_part)
        except NameError:
            result = result_part

mask = (0.09 < result['zspec']) & (result['zspec'] < 0.11)
mass = result[mask]

limits = pyl.logspace(10,15)
sums = pyl.zeros_like(limits)
for i, l in enumerate(limits):
    sums[i] = sum(mass['m200c']>l)

# now we make the volumes
sr = 398 * (pyl.pi/180.)**2
# integrate to the redshift
Vc = quad(astCalc.dVcdz, 0.09, 0.11)[0] * sr

#get hmf data
hmf = pyl.genfromtxt('mVector_PLANCK-SMT .txt')

# plot! -- the 0.7 is because that is H0 used in astLib
pyl.plot(hmf[:,0]/0.7, hmf[:,8]*0.7**3)
pyl.scatter(limits/0.7,sums/Vc)

pyl.loglog()
pyl.xlabel('Mass $(M_\odot)$')
pyl.ylabel('n(>M) Mpc$^{-3}$')
pyl.show()
