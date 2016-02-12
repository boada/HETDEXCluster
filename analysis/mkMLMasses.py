import numpy as np
import h5py as hdf
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import train_test_split
from numpy.lib import recfunctions as rfns

def updateArray(data):
    ''' Adds the results containers to the data product. '''

    newData = np.zeros(data.size)
    data = rfns.append_fields(data, ['ML_pred_1d', 'ML_pred_2d',
        'ML_pred_3d'],
            [newData, newData, newData], dtypes='>f4',
            usemask=False)

    newnewData = np.zeros(data.size, dtype=[('ML_pred_1d_err', '>f4', (2,)),
        ('ML_pred_2d_err', '>f4', (2,)), ('ML_pred_3d_err', '>f4', (2,)),])
    data = rfns.merge_arrays((data, newnewData), usemask=False,
            asrecarray=False, flatten=True)

    return data


### Targeted ###
################
with hdf.File('result_targetedIdeal.hdf5', 'r') as f:
    dset  = f[f.keys()[0]]
    data = dset.value

# add the extra fields
data = updateArray(data)

# You have to clean the data here. This is almost certainly from the fact that
# some of the HALOIDS are repeated at different redshifts. I have a prior on
# the LOSVD calculation which will limit the LOSVD to a maxium. Because the
# clusters are so far apart the LOSVD is super high.
mask = (np.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5)
maskedDataT = data[~mask]
badData = data[mask]
trainT, testT = train_test_split(maskedDataT, test_size=0.5)

### Survey ###
##############
with hdf.File('surveyComplete_noRotations.hdf5', 'r') as f:
    dset  = f[f.keys()[0]]
    data = dset.value

# add the extra fields
data = updateArray(data)

# You have to clean the data here. This is almost certainly from the fact that
# some of the HALOIDS are repeated at different redshifts. I have a prior on
# the LOSVD calculation which will limit the LOSVD to a maxium. Because the
# clusters are so far apart the LOSVD is super high.
mask = (np.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5)
maskedDataS = data[~mask]
badData = data[mask]
trainS, testS = train_test_split(maskedDataS, test_size=0.5)

for train, test in zip([trainT, trainS, testT, testS], [testT, testS, trainT,
    trainS]):
    rf = RandomForestRegressor(n_estimators=1000, min_samples_leaf=1, verbose=1)
    X = np.log10(train['M200c'])

############
#### 1d ####
############
    y = np.column_stack([np.log10(train['LOSVD'])])
    rf.fit(y, X)
    obs = np.column_stack([np.log10(test['LOSVD'])])
    mrf = rf.predict(obs)

    data







#############
#### 2d #####
#############
    y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC']])
    rf.fit(y, X)
    obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC']])
    mrf = rf.predict(obs)

##############
##### 3d #####
##############
    y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC'], train['NGAL']])
    rf.fit(y, X)
    obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC'], test['NGAL']])
    mrf = rf.predict(obs)

