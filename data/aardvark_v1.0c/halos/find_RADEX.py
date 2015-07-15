from glob import glob
import h5py as hdf
from numpy import where

files = glob('*.hdf5')

outFile = open('tiles.txt', 'w')

for f in files:
    print f
    with hdf.File(f, 'r') as f:
        dset = f[f.keys()[0]]
        ra = dset['RA']
        dec = dset['DEC']

        if ra.max() > 300. and ra.min() < 100:
            x = where(ra < 200)
            y = where(ra > 200)
            outFile.writelines('%s %s %s %s %s\n' % (f.keys()[0], ra[x].max(),
                ra[y].min(), dec.max(), dec.min()))
        else:
            outFile.writelines('%s %s %s %s %s\n' % (f.keys()[0], ra.max(),
                ra.min(), dec.max(), dec.min()))

        f.close()
outFile.close()
