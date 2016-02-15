import numpy as np
import h5py as hdf
from sklearn.ensemble import RandomForestRegressor
#from sklearn.cross_validation import train_test_split
from numpy.lib import recfunctions as rfns
from itertools import permutations


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


def splitData(data, test_size=0.3):
    def splitList(alist, wanted_parts=1):
        ''' Breaks a list into a number of parts. If it does not divide evenly
        then the last list wil have an extra element.

        '''

        length = len(alist)
        return [alist[i*length // wanted_parts: (i+1)*length // wanted_parts]\
                for i in range(wanted_parts)]

    np.random.shuffle(data)
    sl = splitList(data, int(1/test_size))

    c = permutations(range(int(1/test_size)))

    prev_i = -1
    for i,j, k in c:
        if i == prev_i:
            continue
        else:
            test = sl[i]
            train = np.append(sl[j], sl[k])
        prev_i = i
        #print test
        #print train

        yield train, test

def addMasses(data, generator):
    ''' This does all of the heavy lifting to get the new masses assigned to
    the right places.

    '''

    i = 0
    for train, test in generator:
        rf = RandomForestRegressor(n_estimators=1000, min_samples_leaf=1,
                verbose=1)
        X = np.log10(train['M200c'])

    ############
    #### 1d ####
    ############
        y = np.column_stack([np.log10(train['LOSVD'])])
        rf.fit(y, X)
        obs = np.column_stack([np.log10(test['LOSVD'])])
        mrf = rf.predict(obs)

        data['ML_pred_1d'][test['IDX']] = mrf

        # errors
        print('Calculating Error')
        err_down, err_up = pred_ints(rf, obs, percentile=68)
        data['ML_pred_1d_err'][test['IDX']] = np.column_stack([err_down,
            err_up])

    #############
    #### 2d #####
    #############
        y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC']])
        rf.fit(y, X)
        obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC']])
        mrf = rf.predict(obs)

        data['ML_pred_2d'][test['IDX']] = mrf
        # errors
        print('Calculating Error, 2d')
        err_down, err_up = pred_ints(rf, obs, percentile=68)
        data['ML_pred_2d_err'][test['IDX']] = np.column_stack([err_down,
            err_up])
    ##############
    ##### 3d #####
    ##############
        y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC'],
            train['NGAL']])
        rf.fit(y, X)
        obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC'],
            test['NGAL']])
        mrf = rf.predict(obs)

        data['ML_pred_3d'][test['IDX']] = mrf
        # errors
        print('Calculating Error, 3d')
        err_down, err_up = pred_ints(rf, obs, percentile=68)
        data['ML_pred_3d_err'][test['IDX']] = np.column_stack([err_down,
            err_up])

        print(i)
        i+=1

    return data

def pred_ints(model, X, percentile=68):
    ''' Calculates the prediction intervals of the estimators. '''

    err_down = []
    err_up = []
    for x in range(len(X)):
        preds = []
        for pred in model.estimators_:
            try:
                preds.append(pred.predict(X[x][:,np.newaxis]))
            except ValueError:
                preds.append(pred.predict(X[x].reshape(1,-1)))
        err_down.append(np.percentile(preds, (100 - percentile) / 2. ))
        err_up.append(np.percentile(preds, 100 - (100 - percentile) / 2.))

    return err_down, err_up


if __name__ == "__main__":

    ### Targeted ###
    ################
    with hdf.File('result_targetedIdeal.hdf5', 'r') as f:
        dset  = f[f.keys()[0]]
        data = dset.value

    # add the extra fields
    data = updateArray(data)

    # You have to clean the data here. This is almost certainly from the fact
    # that some of the HALOIDS are repeated at different redshifts. I have a
    # prior on the LOSVD calculation which will limit the LOSVD to a maxium.
    # Because the clusters are so far apart the LOSVD is super high.

    mask = (np.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5)
    maskedDataT = data[~mask]
    badData = data[mask]

    sl_targeted = splitData(maskedDataT, 0.3)
    data = addMasses(data, sl_targeted)
    with hdf.File('result_targetedIdeal_masses.hdf5', 'w') as f:
        f['predicted masses'] = data
        f.flush()


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

    sl_survey = splitData(maskedDataS, 0.3)
    data = addMasses(data, sl_survey)
    with hdf.File('surveyComplete_noRotations_masses.hdf5', 'w') as f:
        f['predicted masses'] = data
        f.flush()
