import h5py as hdf
import pylab as pyl

def find_indices(bigArr, smallArr):
    from bisect import bisect_left, bisect_right
    ''' Takes the full halo catalog and picks out the HALOIDs that we are
    interested in. Only returns their indexes. It will need to be combined
    with the result catalog in some other way.

    '''

    inds = []
    sortedind = pyl.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    for i, _ in enumerate(smallArr):
        i1 = bisect_left(sortedbigArr, smallArr[i])
        i2 = bisect_right(sortedbigArr, smallArr[i])
        try:
            inds.append(sortedind[i1:i2])
        except IndexError:
            pass
        if i % 100000 ==0:
            print i

    return inds

f = hdf.File('./result_targetedIdeal.hdf5', 'r')
dset = f[f.keys()[0]]
targeted = dset.value
f.close()
f = hdf.File('./redMapper_matched.hdf5', 'r')
dset = f[f.keys()[0]]
RM = dset.value

matched = find_indices(targeted['HALOID'], RM['HALOID'])

cleaned = [[j,i[0]] for j,i in enumerate(matched) if i]
cleaned = pyl.array(cleaned)
#pyl.scatter(pyl.log10(RM['LAMBDA'][cleaned[:,0]]),
#        pyl.log10(targeted['M200c'][cleaned[:,1]]))


f = pyl.figure(1, figsize=(7, 7*(pyl.sqrt(5.)-1.0)/2.0))
ax = f.add_subplot(111)

ax.errorbar(RM['LAMBDA'][cleaned[:,0]], targeted['M200c'][cleaned[:,1]],
        xerr=RM['LAMBDA_err'][cleaned[:,0]], fmt='o', ecolor='0.8')


ax.set_xlabel('Richness')
ax.set_ylabel('Mass $(M_\odot)$')
#pyl.scatter(RM['LAMBDA'][cleaned[:,0]],targeted['M200c'][cleaned[:,1]])


pyl.show()
