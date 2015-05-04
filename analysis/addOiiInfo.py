import pyfits as pyf
import numpy as np
import h5py as hdf
from astLib import vec_astCalc, astStats

with pyf.open('./oii/sdss12_oii_flux_v2.fits') as f:
    sdssData = f[1].data

    # convert to DES magnitudes
    g = sdssData['g'] - 0.104 * (sdssData['g'] - sdssData['r']) + 0.01
    r = sdssData['r'] - 0.102 * (sdssData['g'] - sdssData['r']) + 0.02

    xdat = vec_astCalc.absMag(r, vec_astCalc.dl(sdssData['redshift']))
    ydat = g - r

    bins = [50,50]
    extent = [[-26,-10],[-1,4]]
    _, locx, locy = pyl.histogram2d(xdat, ydat, range=extent, bins=bins)

# now we loop over the catalog galaxies to add the info
with hdf.File('./out1204878_complete.hdf5', 'r+') as f:
    dset = f[f.keys()[0]]
    catg = dset['g']
    catr = dset['r']
    catRedshift = dset['Z']
    catOii = dset['Oii']

    x = vec_astCalc.absMag(catr, vec_astCalc.dl(catRedshift))
    y = catg - catr

    xbin = np.digitize(x, locx)
    ybin = np.digitize(y, locy)

    i = (locx[xbin-1] < xdat) & (xdat < locx[xbin]) & (locy[ybin-1] < ydat) &
        (ydat < locy[ybin])

