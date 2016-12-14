from targetedIdeal_async import *
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky

halo = mkHalo()
truth = mkTruth()
centrals = (truth['CENTRAL'] == True)
maskedTruth = truth[centrals]

f = hdf.File(
    '../data/redMapper/Buzzard-highres_v1.0_redmapper_lgt20_catalog.hdf5', 'r')
dset = f[f.keys()[0]]
RM = dset.value

RMmatch = SkyCoord(ra=RM['RA'] * u.degree, dec=RM['DEC'] * u.degree)
truthMatch = SkyCoord(ra=maskedTruth['RA'] * u.degree,
                      dec=maskedTruth['DEC'] * u.degree)

# this does the catalog matching
idx, d2d, d3d = match_coordinates_sky(RMmatch, truthMatch)

# find the ones that are within 1''
# this is also the index of the matched galaxies in the RM catalog
close = np.where(d2d.arcsec <= 1)
close = close[0]

# these are the indexes of the truth file which have been matched
matchedIndx = idx[close]

# the haloids of the galaxies which have been matched
hids = maskedTruth['HALOID'][matchedIndx]

# these are he indexes of the corresponding halos to the matched galaxies
haloIndx = find_indices(halo['id'], hids)
haloIndx = np.ravel(haloIndx)

# make the results container
results = np.zeros((matchedIndx.size, ),
                   dtype=[('HALOID', '>i8'),
                          ('RA', '>f4'),
                          ('DEC', '>f4'),
                          ('M200c', '>f4'),
                          ('ZSPEC', '>f4'),
                          ('LAMBDA', '>f4'),
                          ('LAMBDA_err', '>f4'), ])

for j, i in enumerate(matchedIndx):
    results['HALOID'][j] = halo['id'][haloIndx[j]]
    results['RA'][j] = maskedTruth['RA'][i]
    results['DEC'][j] = maskedTruth['DEC'][i]
    results['LAMBDA'][j] = RM['LAMBDA_CHISQ'][close[j]]
    results['LAMBDA_err'][j] = RM['LAMBDA_CHISQ_E'][close[j]]
    results['ZSPEC'][j] = halo['zspec'][haloIndx[j]]
    results['M200c'][j] = halo['m200c'][haloIndx[j]] / 0.70

try:
    os.remove('redMapper_matched.hdf5')
except OSError:
    pass
with hdf.File('redMapper_matched.hdf5', 'w') as f:
    f['redMapper_matched'] = results
