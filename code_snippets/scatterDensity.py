import matplotlib.pyplot as plt
import numpy as np
import scipy

#histogram definition
xyrange = [[0, 1], [0, 1]]  # data range
bins = [25, 25]  # number of bins
thresh = 3  # density threshold

#data definition
N = 1e2
xdat, ydat = np.random.normal(size=N), np.random.normal(1, 0.6, size=N)

# histogram the data
hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
posx = np.digitize(xdat, locx)
posy = np.digitize(ydat, locy)

#select points within the histogram
ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
hhsub = hh[posx[ind] - 1,
           posy[ind] - 1]  # values of the histogram where the points are
xdat1 = xdat[ind][hhsub < thresh]  # low density points
ydat1 = ydat[ind][hhsub < thresh]
hh[hh < thresh] = np.nan  # fill the areas with low density by NaNs

plt.scatter(xdat1, ydat1, alpha=0.7)
plt.imshow(hh.T,
           #cmap='jet',
           extent=np.array(xyrange).flatten(),
           interpolation='none')
plt.colorbar()
plt.show()
