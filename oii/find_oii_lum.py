import pylab as pyl
import pyfits as pyf
from astLib import astCalc, astStats

def train(xdat, ydat):
    bins = [50,50]
    thresh = 3

    hh, locx, locy = pyl.histogram2d(xdat, ydat, range=[[-1,4],[12,22]],
            bins=bins)
    posx = pyl.digitize(xdat, locx)
    posy = pyl.digitize(ydat, locy)

    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    # values of histogram with points in the bins.
    hhsub = hh[posx[ind] - 1, posy[ind] - 1]

    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = pyl.nan # fill the areas with low density by NaNs

    return locx, locy, hh.T

f = pyf.open('sdss12_oii_flux.fits')
data = f[1].data
xdat = data['g'] - data['r']
ydat = data['r']

# make the training data
locx, locy, _ = train(xdat, ydat)

# make flux data
dl = [astCalc.dl(z) for z in data['redshift']]
dl = pyl.array(dl)
flux = data['oii_3726_flux']*4*pyl.pi*dl**2

#####
### THIS IS THE PART THAT GETS LOOPED OVER WHEN IT COMES TO THAT ###
#####

x = 0.49785233
y = 17.706039
# find the bins of the color/mag point of interest
xbin = pyl.digitize([x], locx)
ybin = pyl.digitize([y], locy)

# find all of the points inside the bin we are interested ind
i = (locx[xbin-1] < xdat) & (xdat < locx[xbin]) & (locy[ybin-1] < ydat) & (ydat
        < locy[ybin])

px, x = pyl.histogram(pyl.log10(flux[i]), bins=20, normed=True)
x = pyl.linspace(x[0], x[-1], len(px))

s = astStats.slice_sampler(px, N=1, x=x)

print 'log oii lum', s







