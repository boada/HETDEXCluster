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
    halos = find_indices(data['HALOID'], hids)

    # find the number of galaxies present in all of the halos
    l =[]
    for h in halos: l.append(len(h))

    # make some bins
    bins = pyl.arange(10,110,10)

    # put the lengths into the bins
    binned = pyl.digitize(l, bins=bins)

    for i in range(1, bins.shape[0]):
        mask = (binned == i)

        x = pyl.log(data['LOSVD'][mask]/(data['VRMS'][mask]/pyl.sqrt(3)))
        x = x[~pyl.isnan(x)]

        mean = pyl.mean(x)
        variance = pyl.var(x)
        sigma = pyl.sqrt(variance)
        xr = pyl.linspace(-5, 5,100)
        pyl.plot(xr,normpdf(xr,mean,sigma), label=str(i), c=pyl.cm.jet(i/10.))

        #pyl.hist(x, bins=20, normed=True, histtype='step')

    pyl.show()
