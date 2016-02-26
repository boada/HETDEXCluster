import numpy as np
import h5py as hdf
from sklearn.metrics import median_absolute_error, mean_squared_error

def runningStatistic(stat, true, pred):
    bins = np.arange(11.5,16,0.5)
    indx = np.digitize(true, bins) - 1
    binNumber = len(bins)

    running = []
    for k in xrange(binNumber):
        try:
            s = stat(true[indx==k], pred[indx==k])
        except ValueError:
            s = np.nan
        print '%.3f' % (s), '&',
        running.append(s)
    print ''
    return running

### Targeted ###
################
with hdf.File('./result_targetedRealistic_Probmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'Prob_pred_1d', 'Prob_pred_2d', 'Prob_pred_3d']
# filter bad values
mask = (target['Prob_pred_1d'] != 0)
target = target[mask]

### Survey ###
##############
with hdf.File('./surveyCompleteRealistic_Probmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey = dset['M200c', 'MASS', 'Prob_pred_1d', 'Prob_pred_2d', 'Prob_pred_3d']
# filter bad values
mask = (survey['Prob_pred_1d'] != 0)
survey = survey[mask]

for d in [target, survey]:
    print('power law MAE')
    running = runningStatistic(median_absolute_error, np.log10(d['M200c']),
            np.log10(d['MASS']))
    # print('RMSE')
    # running = np.sqrt(runningStatistic(mean_squared_error, np.log10(d['M200c']),
    #         np.log10(d['MASS'])))

############
#### 1d ####
############
    print('1d')
    running = runningStatistic(median_absolute_error, np.log10(d['M200c']),
            d['Prob_pred_1d'])
    # print('RMSE')
    # running = np.sqrt(runningStatistic(mean_squared_error, np.log10(d['M200c']),
    #         d['Prob_pred_1d']))

#############
#### 2d #####
#############
    print('2d')
    running = runningStatistic(median_absolute_error, np.log10(d['M200c']),
            d['Prob_pred_2d'])
    # print('RMSE')
    # running = np.sqrt(runningStatistic(mean_squared_error, np.log10(d['M200c']),
    #         d['Prob_pred_2d']))

##############
##### 3d #####
##############
    print('3d')
    running = runningStatistic(median_absolute_error, np.log10(d['M200c']),
            d['Prob_pred_3d'])
    # print('RMSE')
    # running = np.sqrt(runningStatistic(mean_squared_error, np.log10(d['M200c']),
    #         d['Prob_pred_3d']))

