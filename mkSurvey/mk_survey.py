from astLib.astCoords import shiftRADec
import numpy as np


def mk_ifus(RA, DEC):
    ''' Generates the IFU arrays, and puts them in RA/DEC space. Returns the
    lower left corner for each of the IFUs. The general shape is 78 ifus which
    are in a sort of octagon pattern, with the center 6 removed for some
    reason.

    '''

    # Generate the 10 x 10 grid
    ifus = np.arange(10)
    xgrid, ygrid = np.meshgrid(ifus, ifus)

    # make ifu mask
    mask = [6, 8, 8, 10, 10]

    for idx, i in enumerate(mask + mask[::-1]):
        edge = (10 - i) / 2
        if not edge:
            pass
        else:
            xgrid[idx][:edge] = -1
            xgrid[idx][-edge:] = -1
            ygrid[idx][:edge] = -1
            ygrid[idx][-edge:] = -1

    # mask out the center 6
    ygridt = ygrid.T
    ygridtr = ygridt.ravel()
    ygridtr[43:46] = -1
    ygridtr[53:56] = -1

    # now reshape it back
    ygridt = ygridtr.reshape(10, 10)
    ygrid = ygridt.T

    mask = ygrid != -1

    x = [shiftRADec(RA, DEC, i * 100, 0)[0] for i in xgrid[mask].ravel()]
    y = [shiftRADec(RA, DEC, 0, i * 100)[1] for i in ygrid[mask].ravel()]
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


def mk_pointings(startRA, startDEC, maxRA, maxDEC):
    ''' Makes the pointings for each of the large IFU grids. Tries to make a
    42x10 degree area using 115x27 pointings. Really, they can be anywhere on
    the sky, but they are in a single region at the moment. Returns in the
    lower left corner of each of the pointings.

    '''

    # figure out the number of pointings that we need.
    numRA = 0
    while 1:
        x = shiftRADec(startRA, startDEC, numRA * 984.4, 0)[0]
        if x >= maxRA:
            break
        else:
            numRA += 1
    numDEC = 0
    while 1:
        x = shiftRADec(startRA, startDEC, 0, numDEC * 984.4)[1]
        if x >= maxDEC:
            break
        else:
            numDEC += 1

    print 'number of pointings', numRA, numDEC

    coords = np.asarray(
        [(x, y) for x in xrange(numRA) for y in xrange(numDEC)])
    # 1320'' = 22', the width of the pointings.
    # 984.4'', the width of the IFU grid plus one gap.
    #x = [shiftRADec(startRA, startDEC, i*1320, 0)[0] for i in coords[:,0]]
    #y = [shiftRADec(startRA, startDEC, 0, i*1320)[1] for i in coords[:,1]]
    x = [shiftRADec(startRA, startDEC, i * 984.4, 0)[0] for i in coords[:, 0]]
    y = [shiftRADec(startRA, startDEC, 0, i * 984.4)[1] for i in coords[:, 1]]

    return x, y


def gen_pointings(startRA, startDEC, maxRA, maxDEC):
    ''' Takes the pointings made in mk_pointings and adds the appropriate
    height/width. Gives the lower left and upper right corner of each pointing.
    Is a generate so it can be called in a loop for each survey made.

    '''

    dec = startDEC
    while dec < maxDEC:
        ra = startRA
        while ra < maxRA:
            yield ra, dec
            ra = shiftRADec(ra, dec, 984.4, 0)[0]
        dec = shiftRADec(ra, dec, 0, 984.4)[1]
