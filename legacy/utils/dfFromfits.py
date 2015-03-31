import pandas as pd
from astropy.table import Table

def read_fits_as_dataframe(filename):
# This is the way to read the FITS data into a numpy structured array
# (using astropy.io.fits.getdata didn't work out of the box
# because it gives a FITSRec)
    table = Table.read(filename)
    data = table._data
# Fix byte order.
# See
# https://github.com/astropy/astropy/issues/1156
    data = data.byteswap().newbyteorder()
    try:
        df = pd.DataFrame.from_records(data)
    except Exception:
        df = pd.DataFrame()
        for k in table.keys():
            try:
                df[k] = data[k]
            except ValueError:
                pass
    return df
