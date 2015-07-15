import pylab as pyl
from mk_survey import *
from matplotlib.patches import Rectangle as rec
from matplotlib.patches import RegularPolygon as RP

ra = 0
dec = 0

if ra > 360:
    ra -= 360.

fig = pyl.figure()
ax = fig.add_subplot(111)

for i, p in enumerate(gen_pointings(ra, dec)):
    print i,p

    if i == 7:
        for i in gen_ifus(p[0], p[1]):
            rec1 = rec((i[0], i[1]), i[2]-i[0], i[3]-i[1], color='#188487')
            ax.add_patch(rec1)
    else:
        gi = gen_ifus(p[0], p[1])
        i1 = gi.next()
        x,y = shiftRADec(p[0], p[1], 480., 480.)

        cir1 = RP((x, y), 8, 480*pyl.cos(pyl.pi/8)**-1/3600., zorder=0,
                fc='none', linestyle='dashed', orientation=pyl.pi/8.)
    ax.add_patch(cir1)

pyl.draw()
pyl.show()


