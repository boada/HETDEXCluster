import pylab as pyl
import h5py as hdf
from scipy import stats
from astLib import astStats
from astroML.density_estimation import bayesian_blocks


def error(true, pred, mu):
    ''' Unused, but kept to see how I did it when I wasn't using the Scipy
    functions. Calculates the error on the mean.

    '''
    print(true.size, end=' ')
    if true.size > 1:
        var = pyl.sum((pred - true - mu)**2) / (true.size - 1)
        sem = pyl.sqrt(var / true.size)
        return sem
    elif true.size == 1:
        return 0
    else:
        return pyl.nan


def bias(true, pred):
    ''' unused, but calculates the mean bias. '''

    if true.size > 0:
        return pyl.sum(pred - true) / true.size
        #return pyl.median(true)
    else:
        return pyl.nan


def runningStatistic(stat, true, pred, **kwargs):
    ''' b = bias and s = uncertainty on that bias '''

    bins = pyl.arange(11.5, 16, 0.5)
    indx = pyl.digitize(true, bins) - 1
    binNumber = len(bins)

    runningb = []
    runnings = []
    for k in range(binNumber):
        #print true[indx==k].size,
        b = pyl.mean(pred[indx == k] - true[indx == k])
        s = stats.sem(pred[indx == k] - true[indx == k])
        #print '$%.2f\pm{%.4f}$ &' % (b,s)
        try:
            mean, var, std = stats.mvsdist(pred[indx == k] - true[indx == k])
            print('$%.2f\pm{%.2f}$ &' % (std.mean(), std.std()), end=' ')
        except ValueError:
            print('$%.2f\pm{%.2f}$ &' % (pyl.nan, pyl.nan), end=' ')
        runningb.append(b)
        runnings.append(s)
    print('')
    return runningb, runnings


def calc_err(pred, true):
    return (pred - true) / true


with hdf.File('./targetedRealistic_masses.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    mill = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']

with hdf.File('./targetedRealistic_masses_Buzzard.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    buzz = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']

f = pyl.figure(1, figsize=(7 * (pyl.sqrt(5.) - 1.0) / 2.0, 7))
ax = f.add_subplot(211)
axs = f.add_subplot(212)

for d, c, style, zo in zip([mill, buzz], ['#e24a33', '#467821'], ['-', '--'],
                           [1, 2, 0]):

    ##############
    ##### 3d #####
    ##############
    print('3d')
    bins = bayesian_blocks(pyl.log10(d['M200c']), p0=0.01)
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        d['ML_pred_3d'],
        pyl.percentile,
        binNumber=bins,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    ax.fill_between(y_[0],
                    quants[:, 2],
                    quants[:, 0],
                    facecolor=c,
                    alpha=0.4,
                    edgecolor=c)

    err = calc_err(10**d['ML_pred_3d'], d['M200c'])
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        err,
        pyl.percentile,
        binNumber=bins,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    axs.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    axs.fill_between(y_[0],
                     quants[:, 2],
                     quants[:, 0],
                     facecolor=c,
                     alpha=0.4,
                     edgecolor=c)

    ##############
    ##### 3d #####
    ##############
    #    print('3d for biases')
    #    running = runningBias(bias, pyl.log10(d['M200c']),
    #            d['ML_pred_3d'], bins=bins)

    #    b = pyl.array(running[0])
    #    s = pyl.array(running[1])
    #    axs.plot(bins, b, style, c=c, zorder=zo)
    #    axs.fill_between(bins, b - s, b + s, facecolor=c, alpha=0.4, edgecolor=c)

    ### Add Legend ###
    ##################
line1 = pyl.Line2D([], [], ls='-', color='#e24a33')
line2 = pyl.Line2D([], [], ls='--', color='#467821')
ax.legend((line1, line2), ('Millennium',
                           'Buzzard', ),
          loc='upper left')

ax.text(14.75,
        13.5,
        '$ML_{\sigma, z, Ngal}$',
        fontsize=18,
        horizontalalignment='center', fontweight='medium')

ax.plot([13, 15.5], [13, 15.5], c='k', zorder=0)
ax.set_ylabel('Log $M_{pred}$ ($M_{\odot}$)')
ax.set_xticklabels([])

axs.set_ylabel('$\epsilon$')
axs.set_xlabel('Log $M_{200c}$ ($M_{\odot}$)')
axs.set_ylim(-0.6, 0.6)
axs.set_yticklabels(axs.get_yticks()[:-2])
axs.axhline(0, c='k', zorder=0)

pyl.show()
