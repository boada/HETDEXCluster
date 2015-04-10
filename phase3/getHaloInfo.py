import h5py as hdf
import numpy as np
from halo_handler import *

def split_list(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list will have an extra element.

    '''
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
        for i in range(wanted_parts) ]

# load the result file
f = hdf.File('out1204878.hdf5', 'r')
dset = f[f.keys()[0]]

# load the coordinate information
RAmax = dset['RA'].max()
RAmin = dset['RA'].min()
DECmax = dset['DEC'].max()
DECmin = dset['DEC'].min()

# Get the unique haloids
haloids = np.unique(dset['HALOID'])

# find the tiles that overlap with our area
#haloTiles = findHaloTile(RAmax, DECmax, RAmin, DECmin)

# loads the halo files
#haloCat = mk_haloCatalog(haloTiles)

d = fix_haloTiles()

print 'loading...',
data = mk_haloCatalog(d['name'])
print 'done'

print 'sorting...',
sortedHids = np.argsort(data['HALOID'])
print 'done'

# split the files up into 10 parts
x = split_list(sortedHids, 19)

for idx, p in enumerate(x):
    print idx
    f = hdf.File('halo'+str(idx).zfill(2)+'.hdf5', 'w')
    f[str(idx)] = data[p]
    f.close()





