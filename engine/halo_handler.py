import numpy as np

data_dir = '/home/boada/scratch/halos/'

#data_dir= '/Users/steven/Projects/desCluster/data/halos/'


def find_indices(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i, _ in enumerate(smallArr):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds.append(sortedind[i1:i2])
        except IndexError:
            pass
        if i % 100000 == 0:
            print(i)

    return inds


def find_indices_bool(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = np.zeros(len(bigArr), dtype=bool)
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i, _ in enumerate(smallArr):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds[sortedind[i1:i2][0]] = True
        except IndexError:
            pass

    return inds


def find_indices_multi(bigArr, smallArr, multi):
    from bisect import bisect_left, bisect_right
    from itertools import chain
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for j, sub in enumerate(multi):
        smallArr2 = smallArr[sub]
        tmp = []
        for i, _ in enumerate(smallArr2):
            i1 = bisect_left(sortedbigArr, smallArr2[i])
            i2 = bisect_right(sortedbigArr, smallArr2[i])
            try:
                tmp.append(sortedind[i1:i2])
            except IndexError:
                pass
        if j % 100000 == 0:
            print(j)

        inds.append(np.array(list(chain(*tmp))))
    return inds


def findRADECmaxmin(ra, dec):
    ''' Figures out the max and min RA/DEC for the observed catalog. Use this
    info to figure out which halo tiles we need to look into to get the rest of
    the infomation. Returns RA/DEC max RA/DEC min.

    '''

    if ra.max() > 300. and ra.min() < 100:
        x = np.where(ra < 200)
        y = np.where(ra > 200)
        return ra[x].max(), dec.max(), ra[y].min(), dec.min()
    else:
        return ra.max(), dec.max(), ra.min(), dec.min()
