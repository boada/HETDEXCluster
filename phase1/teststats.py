#from halo_handler import *
from calc_cluster_props import *
#from data_handler import *
#from mk_survey import fix_tiles
from functions import *
import h5py as hdf
import numpy as np
import multiprocessing
from os import environ
#from astLib.astStats import bootstrap
import time
import myglobals

def make_work(truth, halo, inds):
    #for i in range(len(inds)):
    for i in range(25):
        x = np.where(truth['HALOID'] == halo['HALOID'][inds[i]])
        t = truth[x[0]]
        #print len(t), t.dtype
        a = np.where(t['RHALO'] < t['R200'])
        t = t[a[0]]
        t = rfns.append_fields(t, ['CVRMS', 'CR200', 'CM200'],
            [halo['VRMS'][inds[i]]/np.sqrt(3), halo['R200'][inds[i]],
            halo['M200'][inds[i]]], usemask=False)
        yield t

def make_work2(truth, halo, inds):
    x = np.where(truth['HALOID'] == halo['HALOID'][inds])
    t = truth[x[0]]
    #print len(t), t.dtype
    a = np.where(t['RHALO'] < t['R200'])
    t = t[a[0]]
    t = rfns.append_fields(t, ['CVRMS', 'CR200', 'CM200'],
        [halo['VRMS'][inds]/np.sqrt(3), halo['R200'][inds],
        halo['M200'][inds]], usemask=False)
    return t

def mp_handler_progress(inds):
    p = multiprocessing.Pool(20)
    m = multiprocessing.Manager()
    q = m.Queue()

    #inds = inds[:100]

    args = [(i, q) for i in inds]
    result = p.map_async(mp_worker, args)

    while True:
        if result.ready():
            break
        else:
            size = q.qsize()
            print size
            time.sleep(30)

    return result.get()

def mp_handler(truth, halo, inds):

    p = multiprocessing.Pool(20)
    result = p.map(mp_worker, make_work(truth, halo, inds))
    #result = p.map_async(mp_worker, make_work(truth, halo, inds))
    #p.close()
    #p.join()


    #return result.get()
    return result

def mp_worker(args):
#def mp_worker(t):

    inds, q = args

    t = make_work2(myglobals.truth, myglobals.halo, inds)

    lt = len(t)
    id = t['HALOID'][0]
    q.put(id)
    if lt >= 5:
        #print id
        t = updateArray(t)
        t = findClusterCenterRedshift(t)
        t = findLOSV(t)
        try:
            vduc = np.std(t['LOSV'])
            #vduc_e = bootstrap(t['LOSV'], np.std)
            #vduc_e = ((vduc - vduc_e[0]) - (vduc - vduc_e[1]))/2
            vduc_e = 0

            vdc = np.std(t['LOSV'], ddof=1)
            #vdc_e = bootstrap(t['LOSV'], np.std, ddof=1)
            #vdc_e = ((vdc - vdc_e[0]) - (vdc - vdc_e[1]))/2
            vdc_e = 0

            vdg = np.std(t['LOSV'], ddof=1.5)
            #vdg_e = bootstrap(t['LOSV'], np.std, ddof=1.5)
            #vdg_e = ((vdg - vdg_e[0]) - (vdg - vdg_e[1]))/2
            vdg_e = 0

            vdbi = calcVD_big(t['LOSV'])
            #vdbi_e = bootstrap(t['LOSV'], calcVD_big)
            #vdbi_e = ((vdbi - vdbi_e[0]) - (vdbi - vdbi_e[1]))/2
            vdbi_e = 0

            vdgap = calcVD_small(t['LOSV'])
            #vdgap_e = bootstrap(t['LOSV'], calcVD_small)
            #vdgap_e = ((vdgap - vdgap_e[0]) - (vdgap - vdgap_e[1]))/2
            vdgap_e = 0

            vdt = calcVD_test(t['LOSV'])
            #vdt_e = bootstrap(t['LOSV'], calcVD_test)
            #vdt_e = ((vdt - vdt_e[0]) - (vdt - vdt_e[1]))/2
            vdt_e = 0
            return (id,
                    lt,
                    vduc,
                    vduc_e,
                    vdc,
                    vdc_e,
                    vdg,
                    vdg_e,
                    vdbi,
                    vdbi_e,
                    vdgap,
                    vdgap_e,
                    vdt,
                    vdt_e,
                    t['CVRMS'][0],
                    t['CM200'][0],
                    t['CLUSZ'][0])
        except ZeroDivisionError:
            return (id, lt, -1, -1, -1, -1, -1, -1, -1 , -1 , -1, -1, -1, -1, -1,
                -1, -1)

    else:
        return (id, lt, -1, -1, -1, -1, -1, -1, -1 , -1 , -1, -1, -1, -1, -1,
                -1, -1)

def find_overlap():
    ''' This takes the truth list and halo list and figures out which truth
    files go with which halo files. returns a big array with [(HALOID,
    array[TruthIDs])].
    It's not the greatest function ever but it seems to work.

    '''

    halos = fix_tiles('halos.txt')
    truth = fix_tiles('truth.txt')

    final = []
    for IDH, RA1H, RA2H, DEC1H, DEC2H in zip(halos['name'], halos['RAmax'],
            halos['RAmin'], halos['DECmax'], halos['DECmin']):
#    print IDH, RA2H, RA1H, DEC2H, DEC1H

        tiles = []
        for IDT, RA1T, RA2T, DEC1T, DEC2T in zip(truth['name'], truth['RAmax'],
            truth['RAmin'], truth['DECmax'], truth['DECmin']):

            if (RA2H < RA2T < RA1H or RA2H < RA1T < RA1H) and (DEC2H < DEC2T <
                    DEC1H or DEC2H < DEC1T < DEC1H):
                tiles.append(IDT)
#            print IDT, RA2T, RA1T, DEC2T, DEC1T

        final.append((IDH, np.unique(tiles)))

    return final

if __name__ == "__main__":
    base = '/home/boada/scratch'

    data = find_overlap()
    halo = np.random.choice(range(len(data)))

    print data[halo][0]
    print data[halo][1]

    truth = mk_catalog(data[halo][1])
    halo = mk_haloCatalog((data[halo][0],))

    # Only big clusters:
    x = (halo['NGALS'] >= 20) & (halo['Z'] < 0.5)
    halo = halo[x]

    #mask = np.in1d(halo['HALOID'], truth['HALOID'])
    mask = np.in1d(truth['HALOID'], halo['HALOID'])
    truth = truth[mask]

    x = (1e13 < halo['M200']) & (halo['M200'] <= 1e14)
    h1 = halo[x]
    x = (1e14 < halo['M200']) & (halo['M200'] <= 1e15)
    h2 = halo[x]
    x = (1e15 < halo['M200'])
    h3 = halo[x]

    print len(h1), '1e13 - 1e14 Clusters'
    print len(h2), '1e14 - 1e15 Clusters'
    print len(h3), '1e15 or greater Clusters'

    if len(h2) < 10000:
        h = np.concatenate((h2, h3))
    else:
        print 'OMG'

    h1 = np.random.choice(h1, 50000 - len(h))

    halo = np.concatenate((h1, h))

    # find the common haloids we are able to use
    haloids = np.intersect1d(halo['HALOID'], truth['HALOID'])

    print len(haloids), 'Clusters selcted'

    # pick a smaller subset
    #haloids = np.random.choice(haloids, 50000)

    # filter the array to make it smaller and faster to work on.
    # Only takes galaxies with HALOID that we have randomly selected
    truth = truth[np.in1d(truth['HALOID'], haloids)]

    inds = find_indices(halo['HALOID'], haloids)

    myglobals.truth = truth
    myglobals.halo = halo

    print 'do work'
    #result = mp_handler(truth, halo, inds)
    #result = mp_handler_progress(truth, halo, inds)
    result = mp_handler_progress(inds)

    # write it all to a file
    print 'ID number is', environ['LSB_JOBID']
    f = hdf.File('stat_testing'+str(environ['LSB_JOBID'])+'.hdf5', 'w')
    f['dset'] = result
    f.close()
    print 'done'

