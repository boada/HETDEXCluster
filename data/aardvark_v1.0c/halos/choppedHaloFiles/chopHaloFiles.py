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

d = fix_haloTiles()

print 'loading...',
data = mk_haloCatalog(d['name'])
print 'done'

print 'sorting...',
sortedHids = np.argsort(data['HALOID'])
print 'done'

# split the files up into 19 parts
x = split_list(sortedHids, 15)

for idx, p in enumerate(x):
    print idx
    f = hdf.File('halo'+str(idx).zfill(2)+'.hdf5', 'w')
    f[str(idx)] = data[p]
    f.close()

