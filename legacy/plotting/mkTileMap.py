import pylab as pyl
from matplotlib.patches import Rectangle

fig = pyl.figure()
ax = fig.add_subplot(111, projection='mollweide')
#ax = fig.add_subplot(111)

tiles = pyl.genfromtxt('tiles.txt', dtype=None, names=True)
#tiles = pyl.genfromtxt('unrot_tiles.txt', dtype=None, names=True)

# clean up the RAs
org = 0

ramax = pyl.remainder(tiles['RAmax'] + 360 - org, 360)
ramin = pyl.remainder(tiles['RAmin'] + 360 - org, 360)

ind = ramax > 180
ramax[ind] -= 360
ind = ramin > 180
ramin[ind] -= 360

ramax = -ramax
ramin = -ramin


for name, ramax, ramin, decmax, decmin in zip(tiles['name'], ramax, ramin,
        tiles['DECmax'], tiles['DECmin']):


    print name
    rect1 = Rectangle((pyl.radians(ramin), pyl.radians(decmin)),
            pyl.radians(ramax-ramin), pyl.radians(decmax-decmin), fill=None)
    ax.add_patch(rect1)

    x = ramin + (ramax - ramin)/2
    y = decmin + (decmax - decmin)/2

    pyl.text(pyl.radians(x), pyl.radians(y), name.replace('truth',''),
            ha='center')

pyl.show()
