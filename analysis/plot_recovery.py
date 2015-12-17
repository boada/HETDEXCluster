import pylab as pyl
import h5py as hdf
from scatterDensity import scatterDensity

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
        if i % 1000 ==0:
            print i

    return inds

# load the data
f = hdf.File('./result_FullKnowledge.hdf5', 'r')
dset = f[f.keys()[0]]
truth = dset.value
f.close()

f = hdf.File('./result_FullKnowledge_realistic.hdf5', 'r')
dset = f[f.keys()[0]]
target = dset.value
f.close()

f = hdf.File('./surveyComplete.hdf5', 'r')
dset = f[f.keys()[0]]
survey = dset.value

# find the matching HALOIDS
inds = find_indices(truth['HALOID'], target['HALOID'])
Tinds = pyl.ravel(inds)
inds = find_indices(truth['HALOID'], survey['HALOID'])
Sinds = pyl.ravel(inds)
f, ax = pyl.subplots(2,2, figsize=(6,6), squeeze=True)
ax = ax.ravel()

targetGals = target['NGAL'] /truth['NGAL'][Tinds].astype('float')
surveyGals = survey['NGAL'] /truth['NGAL'][Sinds].astype('float')

scatterDensity(ax[0],truth['ZSPEC'][Tinds], targetGals, scale=pyl.log10,
        bins=[40,40])
scatterDensity(ax[2],truth['ZSPEC'][Sinds], surveyGals, scale=pyl.log10,
        bins=[40,40])

# mass plots
scatterDensity(ax[1],pyl.log10(truth['M200c'][Tinds]), targetGals, scale=pyl.log10,
        bins=[40,40])
scatterDensity(ax[3],pyl.log10(truth['M200c'][Sinds]), surveyGals, scale=pyl.log10,
        bins=[40,40])

# adjsut the plots
ax[0].set_ylabel('Recovery Fraction')
ax[2].set_ylabel('Recovery Fraction')
ax[2].set_xlabel('Redshift')
ax[3].set_xlabel('Log Mass')

ax[0].set_xticklabels([])
ax[1].set_yticklabels([])
ax[1].set_xticklabels([])
ax[3].set_yticklabels([])

ax[0].set_xlim(0,0.5)
ax[2].set_xlim(0,0.5)

ax[1].set_xlim(12,15.5)
ax[3].set_xlim(12,15.5)
ax[1].set_xticks([12,13,14,15])
ax[3].set_xticks([12,13,14,15])


pyl.show()
