from mk_survey import *
from data_handler import *
import multiprocessing
import numpy as np
from numpy.lib import recfunctions as rfns
import h5py as hdf
from os import environ

# shared data... This is where the big catalog goes
import myglobals

def mp_worker(ra, dec, ramax, decmax):
    tile = np.unique(find_tile(ra, dec, data=myglobals.data))
    catalog = mk_catalog(tile)
    # These are the individual IFUs, 96 of them
    for ifu in gen_ifus(ra, dec):
        #print ifu
        try:
            result = np.append(result, apply_mask(ifu[0], ifu[1], ifu[2],
                ifu[3], catalog))
        except NameError:
            result = apply_mask(ifu[0], ifu[1], ifu[2], ifu[3], catalog)
    return result
        #return 0

def mp_worker_wrapper(args):
    return mp_worker(*args)

def mp_handler(ra, dec):
    #result = mp_worker_wrapper(gen_pointings(ra, dec).next())
    return rfns.stack_arrays(result, usemask=False)

if __name__ == "__main__":
# survey bounds
    RAmax = 90
    RAmin = 60
    DECmax = -40
    DECmin = -61

    ra = RAmin
    dec = DECmin
    # make the tiles across the boundary friendly. Dont have to do this in the
    # new version but we'll leave it for the time being
    data = fix_tiles()

    ### Here is the parallization ###
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    #result = p.map(mp_worker_wrapper, izip(gen_pointings(ra, dec),
    #    repeat(catalog, 115*27)))
    result = p.map(mp_worker_wrapper, gen_pointings(ra, dec))

    p.close()
    p.join()


    ###################
    ### DO THE WORK ###
    ###################
    print 'do work'
    result = mp_handler(ra, dec)
    print len(result)
    print result.dtype

    # write it all to a file
    print 'ID number is', environ['LSB_JOBID']
    f = hdf.File('out'+str(environ['LSB_JOBID'])+'.hdf5', 'w')
    f['dset'] = result
    f.close()
