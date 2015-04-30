import pylab as pyl
import h5py as hdf

f = hdf.File('./stat_testing281041.hdf5','r')
dset = f[f.keys()[0]]
result = dset.value

# filter out bad values
z = pyl.log10(result[:,1])
z = pyl.where(pyl.isnan(z) != True)
result = result[z[0]]


#bins = pyl.linspace(5, 2000, 50)
bins = pyl.logspace(0.3, 3, 50)
delta = bins[1] - bins[2]
#y = (result[:,2] - result[:,5])/result[:,5]

import itertools
color_cycle = itertools.cycle(pyl.cm.spectral(pyl.linspace(0,1,6)))

for i in range(2, 8):
    c = color_cycle.next()
    y = (result[:,i] - result[:,8])/result[:,8]
    index = pyl.digitize(result[:,1], bins) -1
    avgs = [pyl.mean(y[index == k]) for k in range(len(bins))]
    pyl.plot(bins-delta/2, avgs, lw=2, color=c)
    z = pyl.where(y > 0)
    index = pyl.digitize(result[:,1][z[0]], bins) -1
    avgs = [pyl.mean(y[z[0]][index == k]) for k in range(len(bins))]
    pyl.plot(bins-delta/2, avgs, lw=2, color=c)

    z = pyl.where(y < 0)
    index = pyl.digitize(result[:,1][z[0]], bins) -1
    avgs = [pyl.mean(y[z[0]][index == k]) for k in range(len(bins))]
    pyl.plot(bins-delta/2, avgs, lw=2, color=c)

pyl.show()
