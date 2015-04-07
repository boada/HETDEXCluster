import pylab as pyl
import pyfits as pyf
from astLib import vec_astCalc

bins = [50,50]
thresh = 3

f = pyf.open('sdss12_oii_flux.fits')
data = f[1].data

xdat = data['g'] - data['r']
ydat = data['r']

# convert to absolute magnitudes
ydat = vec_astCalc.absMag(data['r'], vec_astCalc.dl(data['redshift']))

hh, locx, locy = pyl.histogram2d(xdat, ydat, range=[[-1,4],[-26,-10]], bins=bins)
posx = pyl.digitize(xdat, locx)
posy = pyl.digitize(ydat, locy)

# finds the bins which contain points. posx = 0 for points outside "range"
ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
# values of histogram with points in the bins.
hhsub = hh[posx[ind] - 1, posy[ind] - 1]

xdat1 = xdat[ind][hhsub < thresh] # low density points
ydat1 = ydat[ind][hhsub < thresh]
hh[hh < thresh] = 0 # fill the areas with low density by NaNs

pyl.scatter(xdat1, ydat1, s=20, c='0.8')
pyl.imshow(pyl.log10(hh.T), cmap='gray_r',
        extent=pyl.array([[-1,4],[-26,-10]]).flatten(),
        interpolation='nearest')

pyl.show()
