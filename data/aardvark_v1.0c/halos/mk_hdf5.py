from glob import glob
from astropy import table
import h5py as hdf

files = glob('*.fit')

for f in files:
    print f, '.....',
    t =table.Table.read(f, 'r')
    t.write(f.replace('.fit', '.hdf5'),
            path=f.split('_')[1]+f.split('.')[1])
            #path=f.split('_')[2]+f.split('_')[4].split('.')[1])
    print 'done'
