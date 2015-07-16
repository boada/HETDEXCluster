import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from astLib.astCoords import decimal2hms
import h5py as hdf

def rectangle(m,lon1,lat1,lon2,lat2,ec='0.3',shading=None, step=100, ax=None):
    """Draw a projection correct rectangle on the map. RAmax, DECmax, RAmin,
    DECmin is the order of coordinates. Draws the rectangle counterclockwise
    from the RAmax/DECmax. So, across the zero boundary the lower RA value is the
    maximum. That is, across zero, RA=4 > RA=350.

    """

    class pos:
        x = np.array([])
        y = np.array([])

    def line(lons, lats):
        x,y = (lons, lats)
        pos.x = np.hstack((pos.x,x))
        pos.y = np.hstack((pos.y,y))

    # (lon1,lat1) (RA1, DEC1) -> (lon1, lat2) (RA1, DEC2)
    line([lon1]*step, np.linspace(lat1,lat2,step))

    # (lon1, lat2) (RA1, DEC2) -> (lon2, lat2) (RA2, DEC2)
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
        p = Polygon(zip(d[0], d[1]), facecolor=shading, edgecolor=ec,
                        alpha=0.6, linewidth=2)
        if not ax:
            plt.gca().add_patch(p)
        else:
            ax.add_patch(p)

def format(degree):
    hour = decimal2hms(degree, ':').split(':')[0]
    return hour+'$^h$'

def plot_points(m, data, ax=None):
    for i in range(len(data)):
        number = data['name'][i].split('halos')[-1]
        f2 =\
        hdf.File('../../../data/halos/Aardvark_v1.0_halos_r1_rotated.' + number
                + '.hdf5', 'r')
        dset2 = f2[f2.keys()[0]]
        d2 = dset2['RA','DEC','M200', 'Z']
        x= np.where(d2['M200'] >=1e14)
        d2=d2[x[0]]
        x=np.where(np.round(d2['Z']*10)==3)
        d2=d2[x[0]]
        print number, len(d2)
        x,y = m(d2['RA'], d2['DEC'])
        if not ax:
            plt.gca().scatter(x,y, zorder=5)
        else:
            ax.scatter(x,y,zorder=5, s=3)

        del d2, x, y, dset2

def celestialRA(ra):
    value = abs(360.-ra)
    return value

def main():
    from mpl_toolkits.basemap import Basemap
    #m = Basemap(projection='lcc', llcrnrlon=115, llcrnrlat=-90, urcrnrlon=-70,
    #        urcrnrlat=5,lat_0=-50, lon_0=0)
    #m = Basemap(projection='spstere',boundinglat=5,lon_0=0, celestial=True)
    #m = Basemap(projection='spaeqd',boundinglat=5,lon_0=0, celestial=False)
    golden_mean = (np.sqrt(5)-1.0)/2.0 #aesthetic ratio
    fig_width = 10
    fig_height = fig_width*golden_mean
    fig_size = [fig_width,fig_height]

    m = Basemap(projection='stere', width=111319.444*25,
        height=111319.444*40*golden_mean, lat_ts=-45, lat_0=-50,
        lon_0=75)

    #m = Basemap(projection='stere', width=350000.444*30,
    #        height=300000.444*30*golden_mean, lat_ts=-30, lat_0=-45,
    #        lon_0=75)


    f, ax = plt.subplots(1,1)
    ax = [ax]

    data = np.genfromtxt('buzzard_truth.txt', dtype=None, names=True)

    rectangle(m, data['RAmax'].max(), data['DECmax'].max(), data['RAmin'].min(),
            data['DECmin'].min(), shading='#188487', ec='none', ax=ax[0])

#    for ramax, ramin, decmax, decmin in zip(data['RAmax'], data['RAmin'],
#            data['DECmax'], data['DECmin']):
#            print ramax, decmax, ramin, decmin
            #really is ramin, decmax, ramax, decmin
            #rectangle(m, ramax, decmax, ramin, decmin, shading='#188487',
#            rectangle(m, ramax, decmax, ramin, decmin, shading='none',
#                ec='k', ax=ax[0])

    #rectangle(m, 110, 80, 90, 10, shading='r')
    #rectangle(m, 20, 50, 350, 10, shading='r')
    #rectangle(m, 430, -10, 292, -20, shading='#a60628', ec='none', ax=ax[0])
    #rectangle(m, 430, -10, 292, -20, shading='#a60628', ec='none', ax=ax[1])
    #rectangle(m, 91.439438, -2.336169, 59.538841, -9.7056122,
            #shading='#a60628', ec='none', ax=ax[0])

    #plot_points(m, data, ax=ax[1])

    for z,i in enumerate(ax):
        m.drawparallels(np.arange(-80.,10., 10.), labels=[], ax=i)

        if not z:
            m.drawmeridians(np.arange(3,10)* 15, labels=[1,0,0,0],
            labelstyle='+/-', fmt=format, ax=i)
            m.drawmeridians(np.arange(15,23)* 15, labels=[],
            labelstyle='+/-', fmt=format, ax=i)
        if not z:
            m.drawmeridians(np.arange(3,10)* 15, labels=[],
            labelstyle='+/-', fmt=format, ax=i)
            m.drawmeridians(np.arange(15,23)* 15, labels=[0,1,0,0],
            labelstyle='+/-', fmt=format, ax=i)

        m.drawmeridians(np.arange(10,15)* 15, labels=[0,0,1,0],
        labelstyle='+/-', fmt=format, ax=i)
        m.drawmeridians(np.arange(-2,3)* 15, labels=[0,0,0,1],
        labelstyle='+/-', fmt=format, ax=i)

    plt.draw()
    plt.show()
if __name__ == "__main__":
    main()
