from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def rectangle(m,lon1,lat1,lon2,lat2,c='0.3',shading=None, step=100):
    """Draw a projection correct rectangle on the map."""
    class pos:
        x = np.array([])
        y = np.array([])

    def line(lons, lats):
        x,y = (lons, lats)
        #plt.plot(x,y,c=c)
        pos.x = np.hstack((pos.x,x))
        pos.y = np.hstack((pos.y,y))

    # (lon1,lat1) (RA1, DEC1) -> (lon1, lat2) (RA1, DEC2)
    line([lon1]*step, np.linspace(lat1,lat2,step))

    # (lon1, lat2) (RA2, DEC2) -> (lon2, lat2) (RA1, DEC2)
    if abs(lon1) - abs(lon2) < 0:
        line(np.linspace(lon1,0,step), [lat2]*step)
        line(np.linspace(360,lon2,step), [lat2]*step)
    else:
        line(np.linspace(lon1,lon2,step), [lat2]*step)

    # (lon2, lat2) (RA2, DEC2) -> (lon2, lat1) (RA2, DEC1)
    line([lon2]*step, np.linspace(lat2,lat1,step))

    # (lon1, lat1) (RA1, DEC1) -> (lon2, lat1) (RA2, DEC1)
    if abs(lon2) - abs(lon1) > 0:
        line(np.linspace(lon2,360,step), [lat1]*step)
        line(np.linspace(0,lon1,step), [lat1]*step)
    else:
        line(np.linspace(lon2,lon1,step), [lat1]*step)

    if shading:
        d = m(pos.x, pos.y)
        #m.plot(d[0], d[1], lw=2)
        p = Polygon(zip(d[0], d[1]),facecolor=shading, edgecolor='k',
                        alpha=0.5,linewidth=2)
        plt.gca().add_patch(p)


#m = Basemap(projection='moll', lon_0=0, celestial=True)
#m = Basemap(projection='lcc', llcrnrlon=340, llcrnrlat=0, urcrnrlon=40,
#        urcrnrlat=60, lat_0=50, lon_0=0, celestial=True)
#m = Basemap(projection='lcc', llcrnrlon=115, llcrnrlat=-90, urcrnrlon=-70,
#        urcrnrlat=5,lat_0=-50, lon_0=0)
#m = Basemap(llcrnrlon-180,llcrnrlat=5,urcrnrlon=90,urcrnrlat=10,  projection='lcc',lat_0=-50,lon_0=0)
m = Basemap(projection='spstere',boundinglat=5,lon_0=0, celestial=True)


data = np.genfromtxt('tiles.txt', dtype=None, names=True)
for ramax, ramin, decmax, decmin in zip(data['RAmax'], data['RAmin'],
        data['DECmax'], data['DECmin']):

        rectangle(m, ramax, decmax, ramin, decmin, shading='r')

#rectangle(m, 110, 80, 90, 10, shading='r')
#rectangle(m, 20, 50, 350, 10, shading='r')
#rectangle(m, 353, -50, 339, -58, shading='r')

m.drawparallels(np.arange(0.,81.,10.), labels=[0,1,1,0])
m.drawmeridians(np.arange(-180.,181.,20.), labels=[1,0,0,1])
#m.drawparallels(np.arange(-90.,91.,10.), labels=[0,1,1,0])
#m.drawmeridians(np.arange(-71.,115.,10.), labels=[1,0,0,1])
plt.draw()
plt.show()
