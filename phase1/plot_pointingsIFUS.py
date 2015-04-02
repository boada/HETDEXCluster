import pylab as pyl
from mk_survey import *
from matplotlib.patches import Rectangle as rec
from matplotlib.patches import Circle as cir

RAmax = 430
RAmin = 292
DECmax = 0
DECmin = -5
ra = pyl.rand()* ( RAmax - RAmin) + RAmin
dec = pyl.rand()*(DECmax-DECmin)+DECmin
ra = 0
dec = 0

if ra > 360:
    ra -= 360.

fig = pyl.figure()
ax = fig.add_subplot(111)

# draw the ifu boxes
for p in gen_pointings(ra, dec):
    for i in gen_ifus(p[0], p[1]):
        rec1 = rec((i[0], i[1]), i[2]-i[0], i[3]-i[1], color='#a60628')
        ax.add_patch(rec1)

    print p
    gi = gen_ifus(p[0], p[1])
    i1 = gi.next()
    x,y = shiftRADec(i1[0], i1[1], 467., 467.)
    cir1 = cir((x, y), 660.44/3600., zorder=0, fc='#e24a33')
    ax.add_patch(cir1)

pyl.draw()
pyl.show()


