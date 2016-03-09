import numpy as np
import h5py as hdf
from sklearn.metrics import median_absolute_error, mean_squared_error

def scatter2(true, pred, mu):
    if true.size > 0:
        return np.sum((pred - true - mu)**2) /(true.size - 1)
    else:
        return np.nan

def bias(true, pred):
    if true.size > 0:
        return np.sum(pred - true) /true.size
        #return np.median(true)
    else:
        return np.nan
        
def runningStatistic(stat, true, pred, **kwargs):
    bins = np.arange(11.5,16,0.5)
    indx = np.digitize(true, bins) - 1
    binNumber = len(bins)

    runningb = []
    runnings = []
    for k in xrange(binNumber):
        try:
            b = stat(true[indx==k], pred[indx==k], **kwargs)
            s = scatter2(true[indx==k], pred[indx==k], b)
        except ValueError:
            b = np.nan
            s = np.nan
        print '$%.2f\pm{%.2f}$ &' % (b,s),    
        #print '%.3f' % (s), '&',
        runningb.append(b)
        runnings.append(s)
    print ''
    return


### Perfect ###
###############
with hdf.File('./targetedPerfect_MLmasses_realisticOnly.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    perfect = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']
# filter bad values
mask = (perfect['ML_pred_1d'] != 0)
perfect = perfect[mask]

### Targeted ###
################
with hdf.File('./targetedRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']
# filter bad values
mask = (target['ML_pred_1d'] != 0)
target = target[mask]

### Survey ###
##############
with hdf.File('./surveyCompleteRealistic_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']
# filter bad values
mask = (survey['ML_pred_1d'] != 0)
survey = survey[mask]

for d in [perfect, target, survey]:
    print('power law')
    running = runningStatistic(bias, np.log10(d['M200c']),
            np.log10(d['MASS']))

############
#### 1d ####
############
    print('1d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_1d'])

#############
#### 2d #####
#############
    print('2d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_2d'])

##############
##### 3d #####
##############
    print('3d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_3d'])

    print '-----'
