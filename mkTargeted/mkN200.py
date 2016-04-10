from multiprocessing import Pool
import h5py as hdf
import numpy as np
from data_handler import mkN200Data
from halo_handler import find_indices, find_indices_multi
from astLib import astCalc as aca
from astLib import astCoords as aco
import os

# buzzard simulation cosmology
aca.H0 = 70
aca.OMEGA_M0 = 0.286
aca.OMEGA_L0 = 0.714

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

def worker(pos, gals, center, LOSVD, clusz):

    def findSeperationSpatial(data, center, unit='Mpc'):
        ''' Finds the distance to all of the galaxies from the center of the
        cluster in the spatial plane. Returns the values in Mpc.

        '''

        sep = aco.calcAngSepDeg(center[0], center[1], data['RA'],
                data['DEC'])

        if unit == 'Mpc':
            for index, value in enumerate(data['Z']):
                sep[index] *= aca.da(value)/57.2957795131
        elif unit == 'arcsecond':
            data['SEP'] *= 3600

        return sep

    def findR200(LOSVD, avgz):
        ''' This is from Carlberg1997.  Their equation 8.'''
        return np.sqrt(3) * (LOSVD)/(10*aca.H0 * aca.Ez(avgz))

    #print "PID: %d \t Value: %d" % (os.getpid(), pos)
    sep = findSeperationSpatial(gals, center)
    r200c = findR200(LOSVD, clusz)
    # we just want the number not the actual galaxies.
    n200c = np.where(sep < r200c)[0].size

    return pos, r200c, n200c

def cb_func((pos, r200c, n200c)):
    if pos % 1000 == 0:
        print pos
    results['IDX'][pos] = pos
    results['R200c'][pos] = r200c
    results['N200c'][pos] = n200c

if __name__ == "__main__":
    async_worker = AsyncFactory(worker, cb_func)
    halo = mkN200Data('halo')
    truth = mkN200Data('truth')

    # there are no clusters with mass < 2e11 and more than 5 galaxies
    mask = (halo['upid'] == -1) & (halo['m200c'] > 2e11)
    maskedHalo = halo[mask]
    hids, uniqueIdx = np.unique(maskedHalo['id'], return_index=True)

    # now we find all member halos of the most massive halos
    subHalos = find_indices(halo['upid'], hids)

    # now find all the corresponding galaxies
    gals = find_indices_multi(truth['HALOID'], halo['id'], subHalos)

    # how many clusters are we looking at?
    x = [i for i,g in enumerate(gals) if g.size >=5]

    results = np.zeros((len(x),), dtype=[('IDX', '>i4'),
        ('HALOID', '>i8'),
        ('M200c', '>f4'),
        ('R200c', '>f4'),
        ('NGAL', '>i4'),
        ('N200c', '>i4')])

    print('do work', len(x), 'clusters to go!')
    for j,i in enumerate(x):
        center = (maskedHalo['ra'][uniqueIdx[i]],
                maskedHalo['dec'][uniqueIdx[i]])
        losvd = maskedHalo['vrms'][uniqueIdx[i]]/np.sqrt(3)
        clusz = maskedHalo['zspec'][uniqueIdx[i]]

        async_worker.call(j, truth[gals[i]], center, losvd, clusz)
        # update results array
        results['HALOID'][j] = maskedHalo['id'][uniqueIdx[i]]
        results['NGAL'][j] = gals[i].size
        results['M200c'][j] = maskedHalo['m200c'][uniqueIdx[i]]/0.70

    async_worker.wait()

    try:
        os.remove('n200.hdf5')
    except OSError:
        pass
    with hdf.File('n200.hdf5', 'w') as f:
        f['n200'] = results
