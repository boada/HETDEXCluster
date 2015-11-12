import pylab as pyl
import h5py as hdf
from scatterDensity import scatterDensity

# load the data
f = hdf.File('./result_FullKnowledge.hdf5', 'r')
dset = f[f.keys()[0]]
truth = dset.value
f.close()
f = hdf.File('./result_targetedIdeal.hdf5', 'r')
dset = f[f.keys()[0]]
target = dset.value
f.close()
f = hdf.File('./surveyComplete.hdf5', 'r')
dset = f[f.keys()[0]]
survey = dset.value

# only use the results we actually have results for
mask = truth['NGAL'] >=5
maskedTruth = truth[mask]
maskedTarget = target[mask]
maskedSurvey = survey[mask]

f, ax = pyl.subplots(1,2)

targetGals = maskedTarget['NGAL'] /maskedTruth['NGAL'].astype('float')
surveyGals = maskedSurvey['NGAL'] /maskedTruth['NGAL'].astype('float')

scatterDensity(ax[0], pyl.log10(maskedTruth['M200c']), targetGals,
        scale=pyl.log10)
scatterDensity(ax[1], pyl.log10(maskedTruth['M200c']), surveyGals,
        scale=pyl.log10)

pyl.show()
