import numpy as np
import h5py as hdf
from numpy.lib import recfunctions as rfns
from itertools import permutations
from prob_propigator import prob1d
from prob_propigator_3d import prob2d
from prob_propigator_4d import prob3d


def updateArray(data):
    ''' Adds the results containers to the data product. '''

    newData = np.zeros(data.size)
    data = rfns.append_fields(data, ['Prob_pred_1d', 'Prob_pred_2d',
                                     'Prob_pred_3d'],
                              [newData, newData, newData],
                              dtypes='>f4',
                              usemask=False)

    newnewData = np.zeros(data.size,
                          dtype=[('Prob_pred_1d_err', '>f4', (2, )),
                                 ('Prob_pred_2d_err', '>f4', (2, )),
                                 ('Prob_pred_3d_err', '>f4', (2, )), ])
    data = rfns.merge_arrays((data, newnewData),
                             usemask=False,
                             asrecarray=False,
                             flatten=True)

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
    sl = splitList(data, int(1 / test_size))

    c = permutations(range(int(1 / test_size)))

    prev_i = -1
    for i, j, k in c:
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
        print('1d')
        mrf = prob1d(train, test)
        mask = np.where(np.isnan(mrf['MASS']))[0]
        mrf = np.delete(mrf, mask)
        test2 = np.delete(test, mask)
        print(mrf.size, test2.size)

        data['Prob_pred_1d'][test2['IDX']] = mrf['MASS']
        data['Prob_pred_1d_err'][test2['IDX']] = mrf['MASS_err']

        print('2d')
        mrf = prob2d(train, test)
        mask = np.where(np.isnan(mrf['MASS']))[0]
        mrf = np.delete(mrf, mask)
        test2 = np.delete(test, mask)
        print(mrf.size, test2.size)
        data['Prob_pred_2d'][test2['IDX']] = mrf['MASS']
        data['Prob_pred_2d_err'][test2['IDX']] = mrf['MASS_err']

        print('3d')
        mrf = prob3d(train, test)
        mask = np.where(np.isnan(mrf['MASS']))[0]
        mrf = np.delete(mrf, mask)
        test2 = np.delete(test, mask)
        print(mrf.size, test2.size)
        data['Prob_pred_3d'][test2['IDX']] = mrf['MASS']
        data['Prob_pred_3d_err'][test2['IDX']] = mrf['MASS_err']

        print(i)
        i += 1
    return data


if __name__ == "__main__":

    ### Targeted ###
    ################
    with hdf.File('./result_targetedRealistic.hdf5', 'r') as f:
        dset = f[f.keys()[0]]
        data = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
                    'LOSVD_err', 'MASS', 'LOSVD_dist']
    # add the extra fields
    data = updateArray(data)

    # You have to clean the data here. This is almost certainly from the fact
    # that some of the HALOIDS are repeated at different redshifts. I have a
    # prior on the LOSVD calculation which will limit the LOSVD to a maxium.
    # Because the clusters are so far apart the LOSVD is super high.

    mask = ((np.log10(data['LOSVD']) > 3.12) & (data['M200c'] < 10**14.5) |
            (data['LOSVD'] < 50))
    maskedDataT = data[~mask]
    badData = data[mask]

    sl_targeted = splitData(maskedDataT, 0.3)
    data = addMasses(data, sl_targeted)
    # drop the LOSVD_dist field
    data = rfns.drop_fields(data, 'LOSVD_dist')
    with hdf.File('./result_targetedRealistic_Probmasses.hdf5', 'w') as f:
        f['predicted masses'] = data
        f.flush()

    ### Survey ###
    ##############
    with hdf.File('./surveyCompleteRealistic.hdf5', 'r') as f:
        dset = f[f.keys()[0]]
        data = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
                    'LOSVD_err', 'MASS', 'LOSVD_dist']

    # add the extra fields
    data = updateArray(data)

    # You have to clean the data here. This is almost certainly from the fact
    # that some of the HALOIDS are repeated at different redshifts. I have a
    # prior on the LOSVD calculation which will limit the LOSVD to a maxium.
    # Because the clusters are so far apart the LOSVD is super high.

    mask = ((np.log10(data['LOSVD']) > 3.12) & (data['M200c'] < 10**14.5) |
            (data['LOSVD'] < 50))
    maskedDataS = data[~mask]
    badData = data[mask]

    sl_survey = splitData(maskedDataS, 0.3)
    data = addMasses(data, sl_survey)
    # drop the LOSVD_dist field
    data = rfns.drop_fields(data, 'LOSVD_dist')
    with hdf.File('./surveyCompleteRealistic_Probmasses.hdf5', 'w') as f:
        f['predicted masses'] = data
        f.flush()
