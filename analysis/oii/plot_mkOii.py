import pylab as pyl
import pyfits as pyf
from astLib import astCalc
from matplotlib.patches import Rectangle
from itertools import izip
from multiprocessing import Pool, cpu_count

def absMag(mag, dl):
    return astCalc.absMag(mag, dl)

def mp_wrapper(args):
    return absMag(*args)

# start the workers
p = Pool(cpu_count(), maxtasksperchild=10)

with pyf.open('sdss12_oii_flux_v2.fits') as f:
    sdssData = f[1].data

    # convert to DES magnitudes
#    g = sdssData['g'] - 0.104 * (sdssData['g'] - sdssData['r']) + 0.01
#    r = sdssData['r'] - 0.102 * (sdssData['g'] - sdssData['r']) + 0.02

    g = sdssData['g']
    r = sdssData['r']

    dl = pyl.array(p.map(astCalc.dl, sdssData['redshift'], chunksize=200))
    xdat = pyl.array(p.map(mp_wrapper, izip(r, dl), chunksize=200))
    ydat = g - r

    p.close()
    p.join()

    bins = [50,50]
    extent = [[-26,-10],[-1,4]]
    thresh = 3

    _, locx, locy = pyl.histogram2d(xdat, ydat, range=extent, bins=bins)

    # need the Oii luminosity
    lum = sdssData['oii_3726_flux']*4.*pyl.pi*(dl* 3.0857e24)**2. *1e-17

f, ax = pyl.subplots(1,2, figsize=(7, 7*(pyl.sqrt(5.)-1.0)/2.0), squeeze=True)
ax = ax.ravel()

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

# the CMD on the left
ax[0].scatter(xdat1, ydat1, s=10, c='0.8', edgecolor='0.8')
ax[0].imshow(pyl.log10(hh.T), cmap='gray_r',
        extent=pyl.array(extent).flatten(),
        interpolation='nearest')
ax[0].set_xlabel('$M_r$ (mag)')
ax[0].set_ylabel('g-r (mag)')
ax[0].set_xticks([-24, -20, -16, -12])

# add some text lables
ax[0].text(-20, 1.25, '1', color='#467821', fontsize=20)
ax[0].text(-16.75, 0.5, '2', color='#cf4457', fontsize=20)

### add the histograms and little boxes
# the two boxes
xcoord = [-17.706039, -21.017]
ycoord = [0.49785233, 1.24]
colors = ['#cf4457', '#467821']

for x, y, c in zip(xcoord, ycoord, colors):

    # find the bins of the color/mag point of interest
    xbin = pyl.digitize([x], locx)
    ybin = pyl.digitize([y], locy)

    # find all of the points inside the bin we are interested ind
    i = (locx[xbin - 1] < xdat) & (xdat < locx[xbin]) & \
        (locy[ybin - 1] < ydat) & (ydat < locy[ybin])
    ax[1].hist(pyl.log10(lum[i]), bins=20, normed=True, histtype='step', lw=2,
            edgecolor=c)

    # little boxes
    rec = Rectangle((locx[xbin], locy[ybin]), locx[xbin + 1] - locx[xbin],
            locy[ybin + 1] - locy[ybin], lw=2, zorder=10, fc='none', ec=c)
    ax[0].add_patch(rec)

ax[1].set_xlabel('Log $L_{[O_{II}]}$ (erg/s)')
ax[1].set_ylabel('P($L_{[O_{II}]}| M_r,g-r)$')
ax[1].set_xlim(36, 44)
ax[1].text(40, 0.75, '1', color='#467821', fontsize=20)
ax[1].text(38, 0.6, '2', color='#cf4457', fontsize=20)
pyl.show()
