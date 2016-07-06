import pylab as pyl
from matplotlib.patches import Rectangle
import h5py as hdf

### Targeted ###
################
with hdf.File('./result_targetedPerfect.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    #data = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
    #    'LOSVD_err', 'MASS', 'LOSVD_dist']
    data = dset['ZSPEC', 'M200c', 'LOSVD']

mask = ((pyl.log10(data['LOSVD']) > 3.12) & (data['M200c'] < 10**14.5) |
    (data['LOSVD'] < 50))
data = data[~mask]
badData = data[mask]

bins = [25, 25]
extent = [[0.0, 0.5], [pyl.log10(50), 3.12]]
thresh = 3

xdat = data['ZSPEC']
ydat = pyl.log10(data['LOSVD'])

f, ax = pyl.subplots(1, 2, figsize=(7, 7 * (pyl.sqrt(5.) - 1.0) / 2.0),
                     squeeze=True)
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
ax[0].set_xlabel('Redshift')
ax[0].set_ylabel('Log $\sigma$ (km $s^{-1}$)')
#ax[0].set_xticks([-24, -20, -16, -12])

# add some text lables
ax[0].text(0.16, 2.4, '1', color='#467821', fontsize=20)
ax[0].text(0.36, 2.7, '2', color='#cf4457', fontsize=20)

### add the histograms and little boxes
# the two boxes
xcoord = [0.15, 0.35]
ycoord = [2.3, 2.6]
colors = ['#467821', '#cf4457']

for x, y, c in zip(xcoord, ycoord, colors):

    # find the bins of the color/mag point of interest
    xbin = pyl.digitize([x], locx)
    ybin = pyl.digitize([y], locy)

    # find all of the points inside the bin we are interested ind
    i = (locx[xbin - 1] < xdat) & (xdat < locx[xbin]) & \
        (locy[ybin - 1] < ydat) & (ydat < locy[ybin])
    ax[1].hist(pyl.log10(data['M200c'][i]), bins=10, normed=True,
               histtype='step', lw=2, edgecolor=c)

    # little boxes
    rec = Rectangle((locx[xbin], locy[ybin]), locx[xbin + 1] - locx[xbin],
            locy[ybin + 1] - locy[ybin], lw=2, zorder=10, fc='none', ec=c)
    ax[0].add_patch(rec)

ax[1].set_xlabel('Log $M_{200c}$ ($M_{\odot}$)')
ax[1].set_ylabel('P($M_{200c}| z, Log\, \sigma)$')
ax[1].text(12.25, 1, '1', color='#467821', fontsize=20)
ax[1].text(13.75, 1, '2', color='#cf4457', fontsize=20)
pyl.show()
