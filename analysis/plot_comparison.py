import pylab as pyl
import h5py as hdf
from addHaloInfo import find_indices
from matplotlib.mlab import normpdf

#with hdf.File('out1204878_hetdex.hdf5', 'r') as f:
with hdf.File('out1204878_complete.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value

    mask = (data['M200']/0.72 >= 1e13) & (data['Z'] < 0.5)
    data = data[mask]

hids = pyl.unique(data['HALOID'])
halos = pyl.array(find_indices(data['HALOID'], hids))

# find the number of galaxies present in all of the halos
l =[]
for h in halos: l.append(data['M200'][h[0]]/0.72)
#for h in halos: l.append(len(h))

# make some bins
#bins = pyl.arange(10,110,10)
bins = pyl.logspace(13,15,10)

# put the lengths into the bins
binned = pyl.digitize(l, bins=bins)

for i in range(1, bins.shape[0]):
    mask = (binned == i)

    indices = [h[0] for h in halos[mask]]

    x = pyl.log(data['LOSVD'][indices]/(data['VRMS'][indices]/pyl.sqrt(3)))
    x = x[~pyl.isnan(x)]

    mean = pyl.mean(x)
    variance = pyl.var(x)
    sigma = pyl.sqrt(variance)
    xr = pyl.linspace(-1, 1,100)
    #pyl.axvline(mean, c=pyl.cm.jet(i/10.))
    #pyl.plot(xr,normpdf(xr,mean,sigma), label=str(i), c=pyl.cm.jet(i/10.))
    pyl.plot(xr,normpdf(xr,mean,sigma), label=str(len(x)),
            c=pyl.cm.jet(i/10.))

    #pyl.hist(x, bins=20, normed=True, histtype='step')

pyl.show()
