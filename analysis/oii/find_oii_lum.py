import pylab as pyl
import pyfits as pyf
from astLib import astCalc, astStats
from multiprocessing import Pool, cpu_count
from itertools import izip

def absMag(mag, dl):
    return astCalc.absMag(mag, dl)

def mp_wrapper(args):
    return absMag(*args)

# start the workers
p = Pool(cpu_count(), maxtasksperchild=10)

with pyf.open('sdss12_oii_flux_v2.fits') as f:
    sdssData = f[1].data

    # convert to DES magnitudes
    g = sdssData['g'] - 0.104 * (sdssData['g'] - sdssData['r']) + 0.01
    r = sdssData['r'] - 0.102 * (sdssData['g'] - sdssData['r']) + 0.02

    dl = pyl.array(p.map(astCalc.dl, sdssData['redshift'], chunksize=200))
    xdat = pyl.array(p.map(mp_wrapper, izip(r, dl), chunksize=200))
    ydat = g - r

    p.close()
    p.join()

    bins = [50,50]
    extent = [[-26,-10],[-1,4]]
    _, locx, locy = pyl.histogram2d(xdat, ydat, range=extent, bins=bins)

    # need the Oii luminosity
    lum = sdssData['oii_3726_flux']*4.*pyl.pi*(dl* 3.0857e24)**2. *1e-17

#####
### THIS IS THE PART THAT GETS LOOPED OVER WHEN IT COMES TO THAT ###
#####

x = -17.706039
y = 0.49785233
#x = -21.017
#y = 1.24

# find the bins of the color/mag point of interest
xbin = pyl.digitize([x], locx)
ybin = pyl.digitize([y], locy)

# find all of the points inside the bin we are interested ind
i = (locx[xbin-1] < xdat) & (xdat < locx[xbin]) & (locy[ybin-1] < ydat) & (ydat
        < locy[ybin])

px, x = pyl.histogram(pyl.log10(lum[i]), bins=20, normed=True)
x = pyl.linspace(x[0], x[-1], len(px))

s = astStats.slice_sampler(px, N=1, x=x)

print 'log oii lum', s

