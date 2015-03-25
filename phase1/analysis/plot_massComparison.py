import pylab as pyl
from mpl_toolkits.axes_grid1 import axes_grid
import h5py as hdf
from astLib import astCalc as aca


def mass(vd, z, a):
    return (vd/a)**(3.) * 1e15/(aca.H0 * aca.Ez(z)/100)


f = hdf.File('../results/stat_testing346304.hdf5', 'r')
dset = f[f.keys()[0]]
result = dset.value

# filter out bad values
x = pyl.where(result[:,2] != -1)
x = (result[:,2] != -1) & (result[:,1] >= 20)
result = result[x]

# make the masses
m = [mass(vd, z, 1142.35938246) for vd, z in zip(result[:,8], result[:,-1])]
m = pyl.log10(pyl.asarray(m))
result[:,-2] = pyl.log10(result[:,-2])

# make 1 to 1 line
x = pyl.linspace(12,16)

# plot!
ax1 = pyl.subplot2grid((3,4), (0,0), colspan=2, rowspan=3)
ax1.scatter(result[:,-2], m, c='0.8', edgecolor='0.8')
ax1.plot(x, x, c='k')

ax1.axhspan(12,12.5, fc='#348ABD', ec='w', zorder=0)
ax1.axhspan(12.5,13, fc='#7A68A6', ec='w', zorder=0)
ax1.axhspan(13,13.5, fc='#CF4457', ec='w', zorder=0)
ax1.axhspan(13.5,14, fc='#A60628', ec='w', zorder=0)
ax1.axhspan(14,14.5, fc='#188487', ec='w', zorder=0)
ax1.axhspan(14.5,15, fc='#E24A33', ec='w', zorder=0)

ax1.set_xlim(12.5, 15.5)
ax1.set_ylim(11.9, 15.3)



s = (12 < m) & (m < 12.5)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (2,2))
ax2.hist(r, bins=20, normed=True, fc='#348ABD')

s = (12.5 < m) & (m < 13.)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (2,3))
ax2.hist(r, bins=20, normed=True, fc='#7A68A6')

s = (13. < m) & (m < 13.5)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (1,2))
ax2.hist(r, bins=20, normed=True, fc='#CF4457')

s = (13.5 < m) & (m < 14.)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (1,3))
ax2.hist(r, bins=20, normed=True, fc='#A60628')

s = (14. < m) & (m < 14.5)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (0,2))
ax2.hist(r, bins=20, normed=True, fc='#188487')

s = (14.5 < m) & (m < 15.)
r = result[:,-2][s]
ax2 = pyl.subplot2grid((3,4), (0,3))
ax2.hist(r, bins=20, normed=True, fc='#E24A33')


pyl.show()
