import numpy as np
from astLib import astCalc as aca

aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def calcMass(vd, A1D = 1082.9, alpha=0.3361):
    avgz = 0.0
    return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

def calcLOSVD(M, z, A1D = 1082.9, alpha=0.3361):
    return A1D * (M * (aca.H0 * aca.Ez(z)/100.)/1e15)**alpha

def prob2d(train, test):
    # make the joint probability
    data = np.column_stack((np.log10(train['M200c']), np.log10(train['LOSVD']),
        train['ZSPEC']))

    Ngrid = 41, 21, 5
    grid = [np.linspace(data[:,i].min(), data[:,i].max(), g+1) for i,g in
            enumerate(Ngrid)]

    H, (xbins, ybins, zbins) = np.histogramdd(data, bins=grid)

    # The bin order in the histogram now becomes z, y, x!!!
    H = H.T

    # Now we try to recover, and predict
    expected_mass = np.zeros((test.size,), dtype=[('MASS', '>f4'),
        ('MASS_err', '>f4', (2,))])

    for j, (s, z) in enumerate(zip(test['LOSVD_dist'], test['ZSPEC'])):
        # using the distribution of sigmas we make P(s) and ds
        Ps, ds = np.histogram(np.log10(s), bins=50, density=True)
        centers = (ds[1:] + ds[:-1])/2.0

        # now we have to make P(m|s,z) -- from the training sample
        # will need to do this a bunch of times for each sigma in centers
        iy  = np.digitize(centers, ybins)
        iz = np.digitize([z], zbins)[0]
        # now we make P(m) using P(m|s,z) and P(s)ds from above
        Pm = np.zeros(xbins.size -1)
        Psds = Ps*np.diff(ds)
        for i, b in enumerate(iy):
                try:
                    if H[iz-1, b-1, :].sum() == 0.0:
                        Pm_sz = np.zeros(xbins.size -1)
                    else:
                        Pm_sz = H[iz-1, b-1, :]/ H[iz-1, b-1, :].sum()
                except IndexError:
                    Pm_sz = np.zeros(xbins.size -1)
                #print Pm_sz
                Pm += Pm_sz * Psds[i]

        # now we can calculate the expected mass
        centers = (xbins[:-1] + xbins[1:])/2.
        norm = np.sum(Pm * np.diff(xbins))
        M_expect = np.sum(centers * Pm * np.diff(xbins))/norm
        variance = np.sum((centers-M_expect)**2 * Pm* np.diff(xbins))/norm

        expected_mass['MASS'][j] = M_expect
        expected_mass['MASS_err'][j] = [M_expect - np.sqrt(variance), M_expect
                + np.sqrt(variance)]

    #mask = np.where(np.isnan(expected_mass['MASS']))[0]
    #expected_mass = np.delete(expected_mass, mask)
    #test = np.delete(test, mask)

    return expected_mass

if __name__ == "__main__":
    prob2d()
