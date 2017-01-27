import pylab as pyl
import h5py as hdf
from astLib import astStats


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
        if i % 1000 == 0:
            print(i)

    return inds

# load the data
with hdf.File('../result_targetedPerfect.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    truth = dset['HALOID', 'NGAL', 'M200c', 'ZSPEC']

with hdf.File('../result_targetedRealistic.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    target = dset['HALOID', 'NGAL', 'M200c', 'ZSPEC']
    #target = target[target['ZSPEC'] <= 0.2]

with hdf.File('../surveyCompleteRealistic.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    survey = dset['HALOID', 'NGAL', 'M200c', 'ZSPEC']
    #survey = survey[survey['ZSPEC'] <=0.2]

    # find the matching HALOIDS
inds = find_indices(truth['HALOID'], target['HALOID'])
Tinds = pyl.ravel(inds)
inds = find_indices(truth['HALOID'], survey['HALOID'])
Sinds = pyl.ravel(inds)
f, ax = pyl.subplots(
    1, 2, figsize=(7, 7 * (pyl.sqrt(5.) - 1.0) / 2.0),
    squeeze=True)
ax = ax.ravel()

targetGals = target['NGAL'] / truth['NGAL'][Tinds].astype('float')
surveyGals = survey['NGAL'] / truth['NGAL'][Sinds].astype('float')

# redshift plots first
### Targeted ###
y_ = astStats.runningStatistic(truth['ZSPEC'][Tinds], targetGals,
                               pyl.percentile, binNumber=20, q=[16, 50, 84])

quants = pyl.array(y_[1])
ax[0].plot(y_[0], quants[:, 1], c='#7A68A6')
ax[0].fill_between(y_[0],
                   quants[:, 2],
                   quants[:, 0],
                   facecolor='#7A68A6',
                   alpha=0.4,
                   edgecolor='#7A68A6',
                   label='Targeted')

### Survey ###
y_ = astStats.runningStatistic(truth['ZSPEC'][Sinds], surveyGals,
                               pyl.percentile, binNumber=20, q=[16, 50, 84])

quants = pyl.array(y_[1])
ax[0].plot(y_[0], quants[:, 1], '--', c='#188487')
ax[0].fill_between(y_[0],
                   quants[:, 2],
                   quants[:, 0],
                   facecolor='#188487',
                   alpha=0.4,
                   edgecolor='#188487',
                   label='Survey')

# mass plots
### Targeted ###
y_ = astStats.runningStatistic(
    pyl.log10(truth['M200c'][Tinds]),
    targetGals,
    pyl.percentile,
    binNumber=20,
    q=[16, 50, 84])

quants = pyl.array(y_[1])
ax[1].plot(y_[0], quants[:, 1], c='#7A68A6')
ax[1].fill_between(y_[0],
                   quants[:, 2],
                   quants[:, 0],
                   facecolor='#7A68A6',
                   alpha=0.4,
                   edgecolor='#7A68A6')

### Survey ###
y_ = astStats.runningStatistic(
    pyl.log10(truth['M200c'][Sinds]),
    surveyGals,
    pyl.percentile,
    binNumber=20,
    q=[16, 50, 84])

quants = pyl.array(y_[1])
ax[1].plot(y_[0], quants[:, 1], '--', c='#188487')
ax[1].fill_between(y_[0],
                   quants[:, 2],
                   quants[:, 0],
                   facecolor='#188487',
                   alpha=0.4,
                   edgecolor='#188487')

### Add Legend ###
##################
line1 = pyl.Line2D([], [], ls='-', color='#7A68A6')
line2 = pyl.Line2D([], [], ls='--', color='#188487')
ax[0].legend((line1, line2), ('Targeted', 'Survey'), loc='upper right')

# adjsut the plots
ax[0].set_ylabel('Recovery Fraction')
ax[0].set_xticklabels([0., 0.1, 0.2, 0.3, 0.4])
ax[0].set_xlabel('Redshift')
ax[1].set_xlabel('Log $M_{200c}$')

ax[1].set_yticklabels([])

ax[1].set_xlim(0, 1)

ax[1].set_xlim(12, 15.5)
ax[1].set_xticks([12, 13, 14, 15])

pyl.show()
