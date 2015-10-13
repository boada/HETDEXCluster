from multiprocessing import Pool
import h5py as hdf
import numpy as np
from data_handler import mkTruth, apply_mask
from mk_survey import gen_pointings, mk_ifus
import os
from astLib.astCoords import shiftRADec
class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool(maxtasksperchild=10)

    def call(self,*args, **kwargs):
        return self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()

def worker(pos, data, ifus):
    for ifu1, ifu2 in zip(ifus[0], ifus[1]):
        ifu3 = shiftRADec(ifu1, ifu2, 50, 0)[0]
        ifu4 = shiftRADec(ifu1, ifu2, 0, 50)[1]

        result_part = apply_mask(ifu1, ifu2, ifu3, ifu4, data)
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part

    return pos, result

def cb_func((pos, data)):
    #print "PID: %d \t Pos: %d" % (os.getpid(), pos)
    if pos % 500 == 0:
        print pos

if __name__ == "__main__":
    async_worker = AsyncFactory(worker, cb_func)
    truth = mkTruth()

    # brightness limits
    gmask = truth['g'] < 22.
    Oiimask = truth['Oii'] > 3.5
    zmask = truth['Z'] > 0.4
    # z > 0.4 & Oii limit | g < 22 & z < 0.4 | Oii limit & z < 0.4
    mask = (zmask & Oiimask) | (gmask & ~zmask) | (Oiimask & ~zmask)

    truth = truth[mask]

    # survey bounds
    RAmax = 90
    RAmin = 60
    DECmax = -40
    DECmin = -61

    objs = []

    for i, pointing in enumerate(gen_pointings(RAmin, DECmin, maxRA=RAmax,
        maxDEC=DECmax)):

        ifus = mk_ifus(pointing[0], pointing[1])
        objs.append(async_worker.call(i, truth, ifus))

        #if i == 15:
        #    break

    async_worker.wait()

    print('make results')
    # now we have to build the results
    for obj in objs:
        try:
            results = np.append(results, obj.get()[1])
        except NameError:
            results = obj.get()[1]

    with hdf.File('observations'+str(os.environ['LSB_JOBID'])+'.hdf5',
            'w') as f:
        f['observations'] = results
