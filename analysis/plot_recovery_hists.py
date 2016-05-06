import pylab as pyl
import h5py as hdf
from scipy import stats

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
with hdf.File('./result_targetedPerfect.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    truth = dset['HALOID','NGAL', 'M200c', 'ZSPEC']

with hdf.File('./result_targetedRealistic.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['HALOID','NGAL', 'M200c', 'ZSPEC']

with hdf.File('./surveyCompleteRealistic.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey = dset['HALOID','NGAL', 'M200c', 'ZSPEC']

# find the matching HALOIDS
inds = find_indices(truth['HALOID'], target['HALOID'])
Tinds = pyl.ravel(inds)
inds = find_indices(truth['HALOID'], survey['HALOID'])
Sinds = pyl.ravel(inds)

targetGals = target['NGAL'] /truth['NGAL'][Tinds].astype('float')
surveyGals = survey['NGAL'] /truth['NGAL'][Sinds].astype('float')

fig, axes = pyl.subplots(nrows=1, ncols=2, figsize=(7,
    7*(pyl.sqrt(5.)-1.0)/2.0))
ax1 = axes[0]
ax2 = axes[1]

# Targeted First
d = stats.binned_statistic_2d(truth['ZSPEC'][Tinds],
        pyl.log10(truth['M200c'][Tinds]), targetGals, bins=50,
        range=[[0.0,0.5],[12,15.5]])
extent = [d[2][0], d[2][-1], d[1][0], d[1][-1]]
im = ax1.imshow(d[0],  extent=extent, interpolation='nearest',origin='lower',
        vmin=0, vmax=1)

# Survey
d = stats.binned_statistic_2d(truth['ZSPEC'][Sinds],
        pyl.log10(truth['M200c'][Sinds]), surveyGals, bins=50,
        range=[[0.0,0.5],[12,15.5]])
extent = [d[2][0], d[2][-1], d[1][0], d[1][-1]]
ax2.imshow(d[0],  extent=extent, interpolation='nearest',origin='lower',
        vmin=0, vmax=1)

# add Colorbar
cbar_ax = fig.add_axes([0.85, 0.20, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)

# Adjust things
ax1.set_xlim(12, 15.5)
ax2.set_xlim(12, 15.5)
ax1.set_xticks([12,13,14,15])
ax2.set_xticks([12,13,14,15])
ax2.set_yticklabels([])
ax1.set_xlabel('Log $M_{200c}$')
ax2.set_xlabel('Log $M_{200c}$')
cbar.set_ticks([0,0.2,0.4,0.6, 0.8, 1])

ax1.set_ylabel('Redshift')
cbar_ax.set_ylabel('Recovery Fraction')

pyl.tight_layout()
pyl.subplots_adjust(wspace=0.05)
fig.subplots_adjust(right=0.8)


pyl.show()
