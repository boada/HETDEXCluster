import h5py as hdf
import numpy as np
from glob import glob

def split_list(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list will have an extra element.

    '''
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
        for i in range(wanted_parts) ]

files = glob('*.hdf5')

for fname in files:
    print fname
    with hdf.File(fname, 'r') as f:
        dset = f[f.keys()[0]]
        result_part = dset.value
        try:
            result = np.append(result, result_part)
        except NameError:
            result = result_part

print 'sorting...'
sorted = np.lexsort((result['RA'], result['DEC']))
print 'done'

# split the files up into 20 parts
x = split_list(sorted, 20)

for idx, p in enumerate(x):
    print idx
    with hdf.File('truth'+str(idx).zfill(2)+'.hdf5', 'w') as f:
        f[str(idx)] = result[p]
