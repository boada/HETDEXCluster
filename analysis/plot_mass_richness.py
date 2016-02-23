import h5py as hdf
import pylab as pyl
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
        if i % 100000 ==0:
            print i

    return inds


with hdf.File('./result_targetedPerfect_MLmasses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    targeted = dset.value
with hdf.File('./surveyCompletePerfect_MLmasses.hdf5', 'r') as f:
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

f = pyl.figure(1, figsize=(7*(pyl.sqrt(5.)-1.0)/2.0,7))
ax1s = pyl.subplot2grid((3,1), (2,0))
ax1 = pyl.subplot2grid((3,1), (0,0), rowspan=2)

xerr = RM['LAMBDA_err'][cleanedTargeted[:,0]] /\
        (RM['LAMBDA'][cleanedTargeted[:,0]] * pyl.log(10))
yerr = targeted['ML_pred_3d_err'][cleanedTargeted[:,1]][:,1] -\
        targeted['ML_pred_3d'][cleanedTargeted[:,1]]


ax1.errorbar(pyl.log10(RM['LAMBDA'][cleanedTargeted[:,0]]),
        targeted['ML_pred_3d'][cleanedTargeted[:,1]],
        xerr=xerr,
        yerr=yerr,
        fmt='o', ecolor='0.8', mfc='#7a68a6', capsize=0.0, label='Targeted')

# add fits
# DOES NOT INCLUDE ERRORS!!!
x = pyl.log(RM['LAMBDA'][cleanedTargeted[:,0]]/60.)
y = pyl.log(10**targeted['ML_pred_3d'][cleanedTargeted[:,1]]* 0.7/1e14)
fit = astStats.OLSFit(pyl.column_stack([x,y]))
x = pyl.array([10,110])
lny = fit['intercept'] + fit['slope']*pyl.log(x/60.)
y = pyl.exp(lny) * 1e14/0.7
ax1.plot(x,y, ls='-', c='#7a68a6', zorder=0)

# survey points
xerr = RM['LAMBDA_err'][cleanedSurvey[:,0]] /\
        (RM['LAMBDA'][cleanedSurvey[:,0]] * pyl.log(10))
yerr = targeted['ML_pred_3d_err'][cleanedSurvey[:,1]][:,1] -\
        targeted['ML_pred_3d'][cleanedSurvey[:,1]]


ax1.errorbar(pyl.log10(RM['LAMBDA'][cleanedSurvey[:,0]]),
        survey['ML_pred_3d'][cleanedSurvey[:,1]],
        xerr=xerr,
        yerr=yerr,
        fmt='o', ecolor='0.8', mfc='#188487', capsize=0.0, label='Survey')

# add fits
# DOES NOT INCLUDE ERRORS!!!
x = pyl.log(RM['LAMBDA'][cleanedSurvey[:,0]]/60.)
y = pyl.log(10**survey['ML_pred_3d'][cleanedSurvey[:,1]]* 0.7/1e14)
fit = astStats.OLSFit(pyl.column_stack([x,y]))
x = pyl.array([10,110])
lny = fit['intercept'] + fit['slope']*pyl.log(x/60.)
y = pyl.exp(lny) * 1e14/0.7
ax1.plot(x,y, ls='--', c='#188487', zorder=0)

#ax1.set_xlim(10,110)
#ax1s.semilogx()
#ax1s.set_xlim(10,110)
#ax1.set_ylim(5e12, 2e15)

ax1s.set_xlabel('Log Richness, $\lambda$')
ax1.set_ylabel('Log $M_{pred}$ $(M_\odot)$')

# add the Rykoff2012 relation
x = pyl.array([10,110])
lny = 1.48 + 1.06*pyl.log(x/60.)
y = pyl.exp(lny) * 1e14/0.7
ax1.plot(pyl.log10(x),pyl.log10(y), ls='-.', c='k', zorder=0, label='Rykoff2012')


# now we are going to do the scatter panel on the bottom
bins = pyl.arange(20,120,10)
x = RM['LAMBDA'][cleanedSurvey[:,0]]
y = targeted['ML_pred_3d'][cleanedTargeted[:,1]]
index = pyl.digitize(x, bins) -1
running = [pyl.std(y[index==k]) for k in range(10)]
running = pyl.array(running)
tmp = pyl.where(~pyl.isnan(running))
print(pyl.mean(running[tmp[0]]))
ax1s.plot(bins, running, drawstyle='steps-mid', c='#7a68a6')
ax1s.axhline(pyl.mean(running[tmp[0]]), zorder=0, lw=1)

y = survey['ML_pred_3d'][cleanedSurvey[:,1]]
running = [pyl.std(y[index==k]) for k in range(10)]
running = pyl.array(running)
tmp = pyl.where(~pyl.isnan(running))
print(pyl.mean(running[tmp[0]]))
ax1s.plot(bins, running, drawstyle='steps-post', c='#188487', linestyle='--')
ax1s.axhline(pyl.mean(running[tmp[0]]), ls='--', zorder=0, lw=1)


ax1s.set_ylabel('$\sigma_{M|\lambda}$')

ax1.set_xticklabels([])
ax1s.set_yticks([0.1, 0.3, 0.5])


pyl.show()
