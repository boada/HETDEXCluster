import pylab as pyl
from mk_survey import *
from mk_bettermap import rectangle
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='spstere',boundinglat=5,lon_0=0, celestial=True)
data = pyl.genfromtxt('tiles.txt', dtype=None, names=True)
for ramax, ramin, decmax, decmin in zip(data['RAmax'], data['RAmin'],
    data['DECmax'], data['DECmin']):
    rectangle(m, ramax, decmax, ramin, decmin, shading='r')

m.drawparallels(np.arange(-81.,0., 10.))
m.drawmeridians(np.arange(0, 361.,20.), labels=[1,1,1,1])

# survey bounds
RAmax = 430
RAmin = 292
DECmax = -10
DECmin = -20
# here we are making 25 surveys
for ra, dec in zip(pyl.rand(10)*(RAmax-RAmin)+RAmin,
        pyl.rand(10)*(DECmax-DECmin)+DECmin):
    if ra > 360:
        ra -= 360.
    s = pyl.asarray([i for i in gen_pointings(ra, dec)])
    x = pyl.where(s[:,2] <200)
    y = pyl.where(s[:,0] > 200)
    if x[0].any() and y[0].any():
        rectangle(m, s[:,2][x].max(), s[:,3].max(), s[:,0][y].min(),
                s[:,1].min(), shading='g')
    else:
        rectangle(m, s[:,2].max(), s[:,3].max(), s[:,0].min(), s[:,1].min(),
            shading='g')

pyl.draw()
pyl.show()

