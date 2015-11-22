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

def findNeighbors(LOSVD, data):
    inds = (LOSVD - 10 < data['VRMS']) & (data['VRMS'] < LOSVD+10)


f = hdf.File('./result_targetedIdeal.hdf5', 'r')
dset = f[f.keys()[0]]
data = dset.value

numBins = 40
bins = np.linspace(200,1300, numBins)
posx = np.digitize(data['LOSVD'], bins)
means = np.zeros(numBins-1)
for i in range(1, numBins):
    inds = np.where(posx == i)[0]
    if inds.size > 0:
        p, dv = np.histogram(data['VRMS'][inds], bins=20, density=True)
        center = (dv[:20] + dv[1:])/2.
        mean = np.sum(center * p * np.diff(dv))
        mean = np.mean(data['VRMS'][inds])
        means[i-1] = mean
    else:
        means[i-1] = np.nan

mask = (300 < data['LOSVD']) & (data['LOSVD'] < 1300)
obs = data['LOSVD'][mask]
true = data['VRMS'][mask]

center = (bins[:numBins-1] + bins[1:])/2.
pred = np.zeros(obs.size)
for i in range(obs.size):
    dist = abs(obs[i]-center)

    w = 1/dist**1.5
    pred[i] = np.sum(w*means)/np.sum(w)
