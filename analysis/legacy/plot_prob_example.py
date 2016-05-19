import pylab as pyl
import h5py as hdf
import corner


### Targeted ###
################
with hdf.File('./result_targetedPerfect.hdf5', 'r') as f:
    dset  = f[f.keys()[0]]
    #data = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
    #    'LOSVD_err', 'MASS', 'LOSVD_dist']
    data = dset['ZSPEC', 'M200c', 'LOSVD',]

mask = ((pyl.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5) |
    (data['LOSVD'] < 50))
data = data[~mask]
badData = data[mask]

X = pyl.column_stack([data['ZSPEC'], pyl.log10(data['M200c']),
    pyl.log10(data['LOSVD'])])
f = corner.corner(X, labels=['z', 'Log $M_{200c}$', 'Log $\sigma$'], bins=50,
        smooth=True, fill_contours=True)


# add the bars. The axes are row indexed

axes = f.axes

# histogram labels
axes[0].set_ylabel('P(z)')
axes[0].yaxis.set_label_position('right')
axes[4].set_ylabel('P(Log $M_{200c}$)')
axes[4].yaxis.set_label_position('right')
axes[8].set_ylabel('P(Log $\sigma$)')
axes[8].yaxis.set_label_position('right')

# redshift
axes[0].axvspan(0.15,0.19,facecolor='#a60628', alpha=0.25)
axes[3].axvspan(0.15,0.19,facecolor='#a60628', alpha=0.25)
axes[6].axvspan(0.15,0.19,facecolor='#a60628', alpha=0.25)

# LOSVD
axes[6].axhspan(2.55,2.65,facecolor='#a60628', alpha=0.25)
axes[7].axhspan(2.55,2.65,facecolor='#a60628', alpha=0.25)
axes[8].axvspan(2.55,2.65,facecolor='#a60628', alpha=0.25)

# mass
axes[3].axhspan(13.25,13.5,facecolor='#a60628', alpha=0.25)
axes[7].axvspan(13.25,13.5,facecolor='#a60628', alpha=0.25)
axes[4].axvspan(13.25,13.5,facecolor='#a60628', alpha=0.25)


