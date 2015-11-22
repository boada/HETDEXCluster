import h5py as hdf
import numpy as np
from astLib import astCalc as aca
from astLib import astStats as ast

aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def calcMass(LOSVD, z, A1D=1082.9, alpha=0.3361):
    mass = 1e15 * (aca.H0 * aca.Ez(z)/100.) * (LOSVD/A1D)**(1/alpha)
    return mass

f = hdf.File('./result_targetedIdeal.hdf5', 'r')
dset = f[f.keys()[0]]
data = dset.value

mask = (300 < data['LOSVD']) & (data['LOSVD'] < 1300)
data = data[mask]

masses = np.zeros(data.size)
for i, LOSVD in enumerate(data['LOSVD']):
    inds = (LOSVD - 25 < data['VRMS']) & (data['VRMS'] < LOSVD + 25)
    if data['VRMS'][inds].size >= 20:
        p, dv = np.histogram(data['VRMS'][inds], bins=20, density=True)

        # now we resample the LOSVD distribution
        x = np.linspace(dv[0], dv[-1], p.size)
        t2 = ast.slice_sampler(p, N=1000, x=x)
        # now we have a new LOSVD distribution that is similar to the original

        # now we make a new mass distribution
        m = map(calcMass, t2, data['ZSPEC'][i] * np.ones(1000))

        masses[i] = np.mean(m)

    else:
        pass
