import numpy as np
from astLib import astCalc as aca
from astLib.astStats import slice_sampler
import h5py as hdf
from sklearn.cross_validation import train_test_split

aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def calcMass(vd, A1D = 1082.9, alpha=0.3361):
    avgz = 0.0
    return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

def plot_training(mass, sigma, bins=41):
    gridx = np.linspace(mass.min(), mass.max(), bins+1)
    gridy = np.linspace(sigma.min(), sigma.max(), bins+1)
    H, xbins, ybins = np.histogram2d(mass, sigma, bins=[gridx, gridy])

    return H, xbins, ybins

f = hdf.File('result_targetedIdeal.hdf5', 'r')
dset  = f[f.keys()[0]]
data = dset.value
mask = (np.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5)
maskedData = data[~mask]
badData = data[mask]
train, test = train_test_split(maskedData, test_size=0.3)

# make the joint probability
H, xbins, ybins = plot_training(np.log10(train['M200c']),
        np.log10(train['LOSVD']), bins=50)
H = H.T

# Now we try to recover, and predict
expected_mass = np.zeros((test.size,), dtype=[('MASS', '>f4'),
    ('MASS_err', '>f4', (2,))])

for j, s in enumerate(test['LOSVD_dist']):
    # using the distribution of sigmas we make P(s) and ds
    Ps, ds = np.histogram(np.log10(s), bins=50, density=True)
    # and get new new sigmas from the histogram2d
    centers = (ds[1:] + ds[:-1])/2.0

    # now we have to make P(mt|s) -- from the training sample
    # will need to do this a bunch of times for each sigma in centers
    iy  = np.digitize(centers, ybins) #change this line
    # now we make P(mt) using P(mt|s) and P(s)ds from above
    Pmt = np.zeros(xbins.size -1)
    Psds = Ps*np.diff(ds)
    for i in range(iy.size):
        try:
            Pmt_s = H[iy[i]] / H[iy[i]].sum()
        except IndexError:
            Pmt_s = np.zeros(xbins.size -1)
        Pmt += Pmt_s * Psds[i]

    # now we can calculate the expected mass
    centers = (xbins[:-1] + xbins[1:])/2.
    norm = np.sum(Pmt * np.diff(xbins))
    M_expect = np.sum(centers * Pmt * np.diff(xbins))/norm

    variance = np.sum((centers-M_expect)**2 * Pmt* np.diff(xbins))/norm

    expected_mass['MASS'][j] = M_expect
    expected_mass['MASS_err'][j] = [M_expect - np.sqrt(variance), M_expect
            + np.sqrt(variance)]
#    resamples = slice_sampler(Pmt,x=centers, N=1000)
#    expected_mass['MASS_err'][j] = np.percentile(resamples, [16, 84])

