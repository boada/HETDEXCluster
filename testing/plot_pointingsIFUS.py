import pylab as pyl
from matplotlib.patches import Rectangle as rec
from matplotlib.patches import RegularPolygon as RP


def mk_pointings(startRA, startDEC):
    ''' Makes the pointings for each of the large IFU grids. Tries to make a
    42x10 degree area using 115x27 pointings. Really, they can be anywhere on
    the sky, but they are in a single region at the moment. Returns in the
    lower left corner of each of the pointings.
    '''

    numRA = 4
    numDEC = 3
    coords = np.asarray([(x, y) for x in xrange(numRA) for y in xrange(numDEC)])
    # 1320'' = 22', the width of the pointings.
    x = [shiftRADec(startRA, startDEC, i*960, 0)[0] for i in coords[:,0]]
    y = [shiftRADec(startRA, startDEC, 0, i*960)[1] for i in coords[:,1]]

    return x, y

def gen_pointings(startRA, startDEC):
    ''' Takes the pointings made in mk_pointings and adds the appropriate
    height/width. Gives the lower left and upper right corner of each pointing.
    Is a generate so it can be called in a loop for each survey made.
    '''

    coords = mk_pointings(startRA, startDEC)
    for x, y in zip(coords[0], coords[1]):
        # Take the RA / DEC and add the 22' widths.
        # yields RAmin/DECmin -- RAmax/DECmax
        yield x, y, shiftRADec(x, y, 1320, 0)[0], shiftRADec(x, y, 0, 1320)[1]

def mk_ifus(RA, DEC):
    ''' Generates the IFU arrays, and puts them in RA/DEC space. Returns the
    lower left corner for each of the IFUs.
    '''

    # Called by gen_ifus
    # Generate the 10 x 10 array
    #coords = np.asarray([(x, y) for x in xrange(10) for y in xrange(10)])

    # Generate the 10 x 10 grid
    ifus = np.arange(10)
    xgrid, ygrid = np.meshgrid(ifus, ifus)

    # make ifu mask
    xmask = [6, 8, 8, 10, 10]
    ymask = [6, 8, 8, 10, 10]

    for idx, i in enumerate(xmask+xmask[::-1]):
        edge = (10-i)/2
        if not edge:
            pass
        else:
            xgrid[idx][:edge] = -1
            xgrid[idx][-edge:] = -1

    for idx, i in enumerate(ymask+ymask[::-1]):
        edge = (10-i)/2
        if not edge:
            pass
        else:
            ygrid[idx][:edge] = -1
            ygrid[idx][-edge:] = -1

    # mask out the center 6
    ygridt = ygrid.T
    ygridtr = ygridt.ravel()
    ygridtr[43:46] = -1
    ygridtr[53:56] = -1

    # now reshape it back
    ygridt = ygridtr.reshape(10,10)
    ygrid = ygridt.T


    mask = ygrid != -1

    x = [shiftRADec(RA, DEC, i*100, 0)[0] for i in xgrid[mask].ravel()]
    y = [shiftRADec(RA, DEC, 0, i*100)[1] for i in ygrid[mask].ravel()]


    return x, y

def gen_ifus(RA, DEC):
    ''' Generates the actual IFU widths. Calls mk_ifus to create the initial
    grid and then adds the proper height/width of each. This is the function to
    use when making IFU selections. Is a generator so it can be called in a
    loop for each of the pointings made below.
    '''

    coords = mk_ifus(RA, DEC)
    for x, y in zip(coords[0], coords[1]):
        # Take the RA / DEC and add the 50'' widths.
        # yields RAmin/DECmin -- RAmax/DECmax
        yield x, y, shiftRADec(x, y, 50, 0)[0], shiftRADec(x, y, 0, 50)[1]

ra = 0
dec = 0

if ra > 360:
    ra -= 360.

fig = pyl.figure(1,figsize=(5,5*(pyl.sqrt(5.)-1.0)/2.0))
ax = fig.add_subplot(111)

for i, p in enumerate(gen_pointings(ra, dec, 1,1)):
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


ax.set_xlim(-0.05, 1.1)
ax.set_ylim(-0.05,0.85)

ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
ax.set_xticks([])
ax.set_yticks([])

pyl.draw()
pyl.show()


