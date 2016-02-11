import numpy
import healpy
import pylab
from functions import *

def angToPix(nside, lon, lat):
    """
    this function pixelize the stars
    Input (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta = numpy.radians(90. - lat)
    phi = numpy.radians(lon)
    return healpy.ang2pix(nside, theta, phi)


# Now construct the map
nside=256
npix = healpy.nside2npix(nside)

gr_min=0.2
gr_max=0.4
g_min=20.5
g_max=22

data_temp = data_y2_new
data_temp = select_color(data_temp, data_temp['MAG_PSF_G'] - data_temp['MAG_PSF_R'], gr_min, gr_max)
selected_data_in = select_in(data_temp, 'MAG_PSF_G', g_min, g_max)
selected_data_out =  select_out(data_temp, 'MAG_PSF_G', g_min, g_max)

'''
data_temp = data_y2_new
data_temp = select_color(data_temp, data_temp['WAVG_MAG_PSF_G'] - data_temp['WAVG_MAG_PSF_R'], gr_min, gr_max)
selected_data_in = select_in(data_temp, 'WAVG_MAG_PSF_G', g_min, g_max)
selected_data_out =  select_out(data_temp, 'WAVG_MAG_PSF_G', g_min, g_max)
'''
'''
data_temp = data_y1_new
#data_temp['WAVGCALIB_MAG_PSF_G'] = data_y1_new['WAVGCALIB_MAG_PSF_G'] - data_y1_new['XCORR_SFD98_G']
#data_temp['WAVGCALIB_MAG_PSF_R'] = data_y1_new['WAVGCALIB_MAG_PSF_R'] - data_y1_new['XCORR_SFD98_R']
data_temp = select_color(data_temp, data_temp['WAVGCALIB_MAG_PSF_G'] - data_temp['WAVGCALIB_MAG_PSF_R'], gr_min, gr_max)
selected_data_in = select_in(data_temp, 'WAVGCALIB_MAG_PSF_G', g_min, g_max)
'''

npix = healpy.nside2npix(nside)
pix = angToPix(nside, selected_data_in.field('RA'),selected_data_in.field('DEC'))
m = numpy.histogram(pix, bins=numpy.arange(npix + 1))[0].astype(float)
m[m == 0] = healpy.UNSEEN
median_counts = numpy.median(m[m != healpy.UNSEEN])
m_smooth = healpy.smoothing(m, sigma=numpy.radians(0.25))
#m_smooth = m

# Do the plotting
pylab.figure(num='fig1', figsize=(15, 8))
pylab.clf()
healpy.mollview(m_smooth, fig='fig1', min=0.5 * median_counts, max=2. * median_counts, xsize=3200, title='g=%.1f-%.1f'%(g_min, g_max), coord='CG', unit='Count Density', margins=(0, 0.2, 0, 0.1))
#healpy.mollview(m_smooth, fig='fig1', min=0.5 * median_counts, max=2. * median_counts, xsize=3200, title='g=%.1f-%.1f'%(g_min, g_max), unit='Count Density', margins=(0, 0.2, 0, 0.1))
healpy.graticule()
pylab.ylim(-1, -0.2)
#pylab.savefig('fig/y1_cmap_%.1f_%.1f.png'%(g_min, g_max), dpi=150, bbox_inches='tight')
#pylab.savefig('fig/cmap_%.1f_%.1f_radec.png'%(g_min, g_max), dpi=150, bbox_inches='tight')
pylab.show()
