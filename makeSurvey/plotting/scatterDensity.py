bins = [50,50]
thresh = 3

extent = [[xdat.min(), xdat.max()], [ydat.min(), ydat.max()]]

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
