bins = [50,50]
thresh = 3

extent = [[xdat.min(), xdat.max()], [ydat.min(), ydat.max()]]

hh, locx, locy = histogram2d(xdat, ydat, range=extent, bins=bins)
posx = digitize(xdat, locx)
posy = digitize(ydat, locy)

# finds the bins which contain points. posx = 0 for points outside "range"
ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
# values of histogram with points in the bins.
hhsub = hh[posx[ind] - 1, posy[ind] - 1]

xdat1 = xdat[ind][hhsub < thresh] # low density points
ydat1 = ydat[ind][hhsub < thresh]
hh[hh < thresh] = 0 # fill the areas with low density by NaNs

scatter(xdat1, ydat1, s=20, c='0.8')
imshow(log10(hh.T), cmap='gray_r',
        extent=array(extent).flatten(),
        interpolation='nearest')
