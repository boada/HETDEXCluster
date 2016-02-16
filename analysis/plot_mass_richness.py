import h5py as hdf
import pylab as pyl

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
        if i % 100000 ==0:
            print i

    return inds


with hdf.File('./result_targetedIdeal_Complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    targeted = dset.value
with hdf.File('./surveyComplete_noRotations_Complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey = dset.value
with hdf.File('./redMapper_matched.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    RM = dset.value

matchedTargeted = find_indices(targeted['HALOID'], RM['HALOID'])
matchedSurvey = find_indices(survey['HALOID'], RM['HALOID'])

# only use the ones that have matches
cleaned = [[j,i[0]] for j,i in enumerate(matchedTargeted) if i]
cleanedTargeted = pyl.array(cleaned)
cleaned = [[j,i[0]] for j,i in enumerate(matchedSurvey) if i]
cleanedSurvey = pyl.array(cleaned)

f = pyl.figure(1, figsize=(7, 7*(pyl.sqrt(5.)-1.0)/2.0))
ax = f.add_subplot(111)

# targeted points
massErr_up = pyl.absolute(10**targeted['ML_pred_3d'][cleanedTargeted[:,1]] -\
        10**targeted['ML_pred_3d_err'][cleanedTargeted[:,1]][:,1])
massErr_down = pyl.absolute(10**targeted['ML_pred_3d'][cleanedTargeted[:,1]] -\
        10**targeted['ML_pred_3d_err'][cleanedTargeted[:,1]][:,0])

ax.errorbar(RM['LAMBDA'][cleanedTargeted[:,0]],
        10**targeted['ML_pred_3d'][cleanedTargeted[:,1]],
        xerr=RM['LAMBDA_err'][cleanedTargeted[:,0]], yerr=pyl.row_stack([massErr_up,
            massErr_down]), fmt='o', ecolor='0.8', mfc='#7a68a6', capsize=0.0,
        label='Targeted')

# survey points
massErr_up = pyl.absolute(10**survey['ML_pred_3d'][cleanedSurvey[:,1]] -\
        10**targeted['ML_pred_3d_err'][cleanedSurvey[:,1]][:,1])
massErr_down = pyl.absolute(10**survey['ML_pred_3d'][cleanedSurvey[:,1]] -\
        10**targeted['ML_pred_3d_err'][cleanedSurvey[:,1]][:,0])

ax.errorbar(RM['LAMBDA'][cleanedSurvey[:,0]],
        10**survey['ML_pred_3d'][cleanedSurvey[:,1]],
        xerr=RM['LAMBDA_err'][cleanedSurvey[:,0]], yerr=pyl.row_stack([massErr_up,
            massErr_down]), fmt='o', ecolor='0.8', mfc='#188487',
        capsize=0.0, ms=8,
        label='Survey')

pyl.loglog()
ax.set_xlim(10,110)
ax.set_ylim(5e12, 2e15)

ax.set_xlabel('Richness')
ax.set_ylabel('$M_{pred}$ $(M_\odot)$')

pyl.show()
