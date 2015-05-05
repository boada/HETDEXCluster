import pyfits as pyf
import numpy as np
import h5py as hdf
from astLib import vec_astCalc, astStats
from collections import Counter

with pyf.open('./oii/sdss12_oii_flux_v2.fits') as f:
    sdssData = f[1].data

    # convert to DES magnitudes
    g = sdssData['g'] - 0.104 * (sdssData['g'] - sdssData['r']) + 0.01
    r = sdssData['r'] - 0.102 * (sdssData['g'] - sdssData['r']) + 0.02

    xdat = vec_astCalc.absMag(r, vec_astCalc.dl(sdssData['redshift']))
    ydat = g - r

    bins = [50,50]
    extent = [[-26,-10],[-1,4]]
    _, locx, locy = np.histogram2d(xdat, ydat, range=extent, bins=bins)

# now we loop over the catalog galaxies to add the info
with hdf.File('./out1204878_complete.hdf5', 'r+') as f:
    dset = f[f.keys()[0]]
    catg = dset['g']
    catr = dset['r']
    catRedshift = dset['Z']
    catOii = dset['Oii']

    # need the Oii luminosity
    lum = sdssData['oii_3726_flux']\
        *4*np.pi*vec_astCalc.dl(sdssData['redshift'])**2

    x = vec_astCalc.absMag(catr, vec_astCalc.dl(catRedshift))
    y = catg - catr

    xbin = np.digitize(x, locx)
    ybin = np.digitize(y, locy)

    # find the unique pairs of bins
    pairs = [(x_,y_) for x_,y_ in zip(xbin, ybin)]
    c = Counter(pairs)

    for bins, number in c.items():
        # need this bit to handle the situation of data outside the bin range
        if bins[0] == len(locx) or bins[1] == len(locy):
            lumes = []
        else:
            # find all of the points inside the bin we are interested in
            i = (locx[bins[0]-1] < xdat) & (xdat < locx[bins[0]]) &\
                (locy[bins[1]-1] < ydat) & (ydat < locy[bins[1]])
            lumes = lum[i]

        if len(lumes) > 10:
            # make a histogram of the flux values
            px, x = np.histogram(np.log10(lumes), bins=20, normed=True)

            # resample the distribution
            x = np.linspace(x[0], x[-1], len(px))
            s = astStats.slice_sampler(px, N=number, x=x)

        elif len(lumes) > 1:
            s = np.mean(np.log10(lumes))
        else:
            s = 0.

        # now we have to find all of the bins.
        binMask = (xbin == bins[0]) & (ybin == bins[1])

        catOii[binMask] = s



