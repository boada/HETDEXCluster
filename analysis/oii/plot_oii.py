import pylab as pyl
import pyfits as pyf
from astLib import vec_astCalc

bins = [50, 50]
thresh = 3
extent = [[-26, -10], [-1, 4]]

with pyf.open('./sdss12_oii_flux_v2.fits') as f:
    data = f[1].data

ydat = data['g'] - data['r']
xdat = data['r']

# convert to des magnitudes -- buzzard is in sdss mags
#g = data['g'] - 0.104*(data['g'] - data['r']) - 0.01
#r = data['r'] - 0.102*(data['g'] - data['r']) - 0.02
#ydat = g-r
#xdat = r

# convert to absolute magnitudes
xdat = vec_astCalc.absMag(xdat, vec_astCalc.dl(data['redshift']))

hh, locx, locy = pyl.histogram2d(xdat, ydat, range=extent, bins=bins)
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
        extent=pyl.array(extent).flatten(),
        interpolation='nearest')

pyl.show()
