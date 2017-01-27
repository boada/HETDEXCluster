import numpy as np
import h5py as hdf
from scipy import stats


def error(true, pred, mu):
    ''' Unused, but kept to see how I did it when I wasn't using the Scipy
    functions. Calculates the error on the mean.

    '''
    print(true.size, end=' ')
    if true.size > 1:
        var = np.sum((pred - true - mu)**2) / (true.size - 1)
        sem = np.sqrt(var / true.size)
        return sem
    elif true.size == 1:
        return 0
    else:
        return np.nan


def bias(true, pred):
    ''' unused, but calculates the mean bias. '''

    if true.size > 0:
        return np.sum(pred - true) / true.size
        #return np.median(true)
    else:
        return np.nan


def runningStatistic(stat, true, pred, **kwargs):
    ''' b = bias and s = uncertainty on that bias '''

    bins = np.arange(12.5, 16, 0.5)
    indx = np.digitize(true, bins)
    binNumber = len(bins)

    runningb = []
    runnings = []
    for k in range(1, binNumber):
        #print true[indx==k].size,
        b = np.mean(pred[indx == k] - true[indx == k])
        s = stats.sem(pred[indx == k] - true[indx == k])
        #print '$%.2f\pm{%.4f}$ &' % (b,s)
        try:
            mean, var, std = stats.mvsdist(pred[indx == k] - true[indx == k])
            print('$%.2f\pm{%.3f}$ &' % (std.mean(), std.std()), end=' ')
        except ValueError:
            print('$%.2f\pm{%.3f}$ &' % (np.nan, np.nan), end=' ')
        runningb.append(b)
        runnings.append(s)
    print('')
    return

### Perfect ###
###############
with hdf.File('./targetedPerfect_MLmasses_realisticOnly.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    perfect = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_2d2',
                   'ML_pred_3d']
# filter bad values
mask = (perfect['ML_pred_1d'] != 0)
perfect = perfect[mask]

### Targeted ###
################
with hdf.File('./targetedRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_2d2',
                  'ML_pred_3d']
# filter bad values
mask = (target['ML_pred_1d'] != 0)
target = target[mask]

### Survey ###
##############
with hdf.File('./surveyCompleteRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[list(f.keys())[0]]
    survey = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_2d2',
                  'ML_pred_3d']
# filter bad values
mask = (survey['ML_pred_1d'] != 0)
survey = survey[mask]

for d in [perfect, target, survey]:

    ### Full survey ###
    mean, var, std = stats.mvsdist(np.log10(d['MASS']) - np.log10(d['M200c']))
    s = stats.sem(np.log10(d['MASS']) - np.log10(d['M200c']))
    #print '$%.2f\pm{%.3f}$' % (mean.mean(),s)
    print('$%.2f\pm{%.3f}$' % (std.mean(), std.std()))

    print('power law')
    running = runningStatistic(bias, np.log10(d['M200c']), np.log10(d['MASS']))
    ############
    #### 1d ####
    ############
    print('1d')
    running = runningStatistic(bias, np.log10(d['M200c']), d['ML_pred_1d'])

    #############
    #### 2d #####
    #############
    print('2d')
    running = runningStatistic(bias, np.log10(d['M200c']), d['ML_pred_2d'])

    #############
    #### 2d2 #####
    #############
    print('2d2')
    running = runningStatistic(bias, np.log10(d['M200c']), d['ML_pred_2d2'])

    ##############
    ##### 3d #####
    ##############
    print('3d')
    running = runningStatistic(bias, np.log10(d['M200c']), d['ML_pred_3d'])

    print('-----')
