from glob import glob
import h5py as hdf
from numpy import where

files = glob('*.hdf5')

outFile = open('HidsMinMax.txt', 'w')

for f in files:
    print f
    try:
        f = hdf.File(f, 'r')
        dset = f[f.keys()[0]]
        hid = dset['HALOID']
        outFile.writelines('%s %s %s\n' % (f.keys()[0], hid.min(),
                hid.max())

        f.close()
    except:
        'Something is wrong!'
outFile.close()
