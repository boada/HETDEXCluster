from mk_survey import *
from data_handler import *
import multiprocessing
import numpy as np

def mp_worker((ra, dec, ramax, decmax)):
    tile = np.unique(find_tile(ra, dec, data=data))
    # These are the individual IFUs, 96 of them
    for ifu in gen_ifus(ra, dec):
        try:
            result = np.append(result, apply_mask(ifu[0], ifu[1], ifu[2],
                ifu[3], tile))
        except NameError:
            result = apply_mask(ifu[0], ifu[1], ifu[2], ifu[3], tile)
    return result
        #return 0


def mp_handler(ra, dec):
    ### Here is the parallization ###
    p = multiprocessing.Pool(2)
    result = p.map(mp_worker, gen_pointings(ra, dec))
    p.close()
    p.join()
    return result

if __name__ == "__main__":
# survey bounds
    RAmax = 430
    RAmin = 292
    DECmax = -10
    DECmin = -20

    # here we are making 25 surveys
    # for ra, dec in zip(pyl.rand(1)*(RAmax-RAmin)+RAmin,
    #                pyl.rand(1)*(DECmax-DECmin)+DECmin):
    for ra, dec in zip([346], [-10]):
        if ra > 360:
            ra -= 360.
        # make the tiles across the boundary friendly.
        data = fix_tiles()
        # now we make the individual pointings, 115x27 of them
        for p in gen_pointings(ra, dec):
            try:
                tile = np.append(tile, find_tile(p[0], p[1], data=data))
                tile = np.append(tile, find_tile(p[2], p[3], data=data))
            except NameError:
                tile = find_tile(p[0], p[1], data=data)
                tile = np.append(tile, find_tile(p[2], p[3], data=data))

	tile = np.unique(tile)
        # load the data -- only the tiles we are going to use
#        global catalog
#        catalog = load_tiles(tile)

        # trim down the the data list to make the finding faster
        # only returns the tiles found in the previous step.
        data = np.asarray([d1 for d1 in data if d1['name'] in tile],
                dtype=data.dtype)

        ###################
        ### DO THE WORK ###
        ###################
        result = mp_handler(ra, dec)
        print len(result)
        #result = mp_worker(gen_pointings(ra,dec).next())

