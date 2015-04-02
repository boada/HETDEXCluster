from astLib.astCoords import shiftRADec
import numpy as np

def mk_ifus(RA, DEC):
    ''' Generates the IFU arrays, and puts them in RA/DEC space. Returns the
    lower left corner for each of the IFUs.

    '''

    # Called by gen_ifus
    # Generate the 10 x 10 array
    #coords = np.asarray([(x, y) for x in xrange(10) for y in xrange(10)])

    # Generate the 10 x 11 array
    xgrid, ygrid = np.mgrid[:10, :11]

    # calculate the distance from the center to each point
    cir = (xgrid - 4.5)**2 + (ygrid - 4)**2

    # draw the hexagon with the center 6 points removed. 72 total points
    c = (cir <= 5**2) & (cir >=2)


    # Remove the center 4 boxes
    #coords2 =  np.append(coords[:44], coords[46:54], axis=0)
    #coords2 = np.append(coords2, coords[56:], axis=0)

    # Makes the RA/DEC grid, ***lower left corner***
    #x = [shiftRADec(RA, DEC, i*98.4, 0)[0] for i in coords2[:,0]]
    #y = [shiftRADec(RA, DEC, 0, i*98.4)[1] for i in coords2[:,1]]

    x = [shiftRADec(RA, DEC, i*98.4, 0)[0] for i in xgrid[c]]
    y = [shiftRADec(RA, DEC, 0, i*98.4)[1] for i in ygrid[c]]
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

def mk_pointings(startRA, startDEC):
    ''' Makes the pointings for each of the large IFU grids. Tries to make a
    42x10 degree area using 115x27 pointings. Really, they can be anywhere on
    the sky, but they are in a single region at the moment. Returns in the
    lower left corner of each of the pointings.

    '''

    # make the area 42 x 10 degrees or 115 x 27 pointings
    numRA = 115
    numDEC = 27
    #numRA = 2
    #numDEC = 2
    coords = np.asarray([(x, y) for x in xrange(numRA) for y in xrange(numDEC)])
    # 1320'' = 22', the width of the pointings.
    x = [shiftRADec(startRA, startDEC, i*1320, 0)[0] for i in coords[:,0]]
    y = [shiftRADec(startRA, startDEC, 0, i*1320)[1] for i in coords[:,1]]

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

def find_tile(RA, DEC, data=False):
    ''' Returns the name of the tile that the current pointing is located
    inside of. Should be called twice for each pointing. Once for the lower
    left and once for the upper right, but only if the lower left passes.
    Throws and error if the pointing sits outside the bounds of any tile.

    '''

    if not len(data):
        data = fix_tiles()
        #data = np.genfromtxt('tiles.txt', names=True, dtype=None)
    else:
        pass
    # Find the RA/DEC of the tile
    tileDEC = (data['DECmax'] > DEC) & (DEC > data['DECmin'])
    tileRA = (data['RAmax'] > RA) & (RA > data['RAmin'])
    tile = np.intersect1d(data['name'][tileRA], data['name'][tileDEC])

    if len(tile):
        return tile
    else:
        raise ValueError('Out of RA/DEC bounds!')

def fix_tiles():
    ''' fixes the tile catalog for the tiles that overlap zero. Now the tiles
    will go from RAmin -> 360 and 0 -> RAmax. Still should only return one tile
    if only one tile is covering the data point.

    '''

    data = np.genfromtxt('tiles.txt', names=True, dtype=None)
    # Find the tiles that overlap zero.
    x = np.where(data['RAmax'] < data['RAmin'])
    d2 = data[x]
    d2['RAmax'] = 360.
    d3 = data[x]
    d3['RAmin'] = 0.0
    data = np.delete(data, x)
    data = np.append(data, d2)
    data = np.append(data, d3)

    return data
