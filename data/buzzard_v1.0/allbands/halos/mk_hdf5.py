from glob import glob
from astropy import table
from astropy.io import fits
from os.path import exists

files = glob('*.fit')

for f in files:
    print f, '.....',
    if not exists(f.replace('fit', 'hdf5')):
        try:
            t = table.Table.read(f, 'r')
            t.write(f.replace('.fit', '.hdf5'),
                    path=f.split('_')[1]+'_'+f.split('.')[1])
            print 'Done'
        except IOError:
            data = fits.getdata(f, ignore_missing_end=True)
            t = table.Table(data)
            t.write(f.replace('.fit', '.hdf5'),
                    path=f.split('_')[1]+'_'+f.split('.')[1])
            print 'Done'
    else:
        print 'Already Done'
