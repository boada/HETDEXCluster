import pyfits as pyf
import numpy as np
import h5py as hdf
from astLib import astCalc, astStats
from collections import Counter
from multiprocessing import Pool, cpu_count
from itertools import izip
from numpy.lib import recfunctions as rfns
import os


def absMag(mag, dl):
    return astCalc.absMag(mag, dl)


def mp_wrapper(args):
    return absMag(*args)


def updateArray(data):
    ''' Updates the data array with the new column we are going to be writing
    to. In this specific case it adds an Oii column.

    '''

    print('updating array...')
    newData = -np.ones(data.shape[0])
    data = rfns.append_fields(data, ['Oii', 'g', 'r', 'Z'], [newData, newData,
                                                             newData, newData],
                              dtypes='>f4',
                              usemask=False)

    return data

# start the workers
p = Pool(cpu_count())

with pyf.open('../analysis/oii/sdss12_oii_flux_v2.fits') as f:
    sdssData = f[1].data

    # convert to DES magnitudes
    #g = sdssData['g'] - 0.104 * (sdssData['g'] - sdssData['r']) + 0.01
    #r = sdssData['r'] - 0.102 * (sdssData['g'] - sdssData['r']) + 0.02

    g = sdssData['g']
    r = sdssData['r']

    dl = np.array(map(astCalc.dl, sdssData['redshift']))
    xdat = np.array(map(mp_wrapper, izip(r, dl)))
    ydat = g - r

    bins = [50, 50]
    extent = [[-26, -10], [-1, 4]]
    _, locx, locy = np.histogram2d(xdat, ydat, range=extent, bins=bins)

    # need the Oii luminosity
    lum = sdssData['oii_3726_flux'] * 4. * np.pi * (dl * 3.0857e24)**2. * 1e-17

# now we loop over the catalog galaxies to add the info
SOL = 2.99e5  # speed of light in km/s
with hdf.File('./galaxies.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    velx = dset['velx']
    cz = dset['redshift']
    catRedshift = velx / SOL * (1 + cz) + cz
    catOii = -np.ones(dset.shape[0])

    dl = np.array(map(astCalc.dl, catRedshift))
    catg = dset['g_sdss'] + 5 * np.log10(dl * 1e6) - 5  # g band apparent mag
    catr = dset['r_sdss'] + 5 * np.log10(dl * 1e6) - 5  # r band apparent mag

    x = dset['r_sdss']  # absolute mag
    y = catg - catr

    xbin = np.digitize(x, locx)
    ybin = np.digitize(y, locy)

    # find the unique pairs of bins
    pairs = [(x_, y_) for x_, y_ in zip(xbin, ybin)]
    c = Counter(pairs)

    for bins, number in c.items():
        # need this bit to handle the situation of data outside the bin range
        if bins[0] == len(locx) or bins[1] == len(locy):
            lumes = []
        else:
            # find all of the points inside the bin we are interested in
            mask = (locx[bins[0]-1] < xdat) & (xdat < locx[bins[0]]) &\
                (locy[bins[1]-1] < ydat) & (ydat < locy[bins[1]])
            lumes = lum[mask]

        # now we have to find all of the bins.
        binMask = (xbin == bins[0]) & (ybin == bins[1])

        if len(lumes) > 10:
            # make a histogram of the flux values
            px, x = np.histogram(np.log10(lumes), bins=20, normed=True)

            # resample the distribution
            x = np.linspace(x[0], x[-1], len(px))
            s = astStats.slice_sampler(px, N=number, x=x)

            catOii[binMask] = (10**s)/(4*np.pi*(3.0857e24*dl[binMask])**2 \
                *1e-17)

        elif len(lumes) >= 1:
            s = np.mean(np.log10(lumes))

            catOii[binMask] = (10**s)/(4*np.pi*(3.0857e24*dl[binMask])**2 \
                *1e-17)
        else:
            catOii[binMask] = 0.

    try:
        os.remove('./galaxies_Oii.hdf5')
    except OSError:
        pass

    with hdf.File('./galaxies_Oii.hdf5', 'w') as f2:
        data = updateArray(dset.value)
        data['Oii'] = catOii
        data['g'] = catg
        data['r'] = catr
        data['Z'] = catRedshift
        f2['./galaxies_Oii.hdf5'] = data
        f2.flush()

p.close()
p.join()

#    dset['Oii'] = catOii
#    f.flush()
