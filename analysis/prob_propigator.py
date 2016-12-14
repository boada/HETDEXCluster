import numpy as np
from astLib import astCalc as aca

aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77


def calcMass(vd, A1D=1082.9, alpha=0.3361):
    avgz = 0.0
    return 1e15 / (aca.H0 * aca.Ez(avgz) / 100.) * (vd / A1D)**(1 / alpha)


def plot_training(mass, sigma, bins=41):
    gridx = np.linspace(mass.min(), mass.max(), bins + 1)
    gridy = np.linspace(sigma.min(), sigma.max(), bins + 1)
    H, xbins, ybins = np.histogram2d(mass, sigma, bins=[gridx, gridy])

    return H, xbins, ybins


def prob1d(train, test):

    # make the joint probability
    data = np.column_stack((np.log10(train['M200c']),
                            np.log10(train['LOSVD'])))

    Ngrid = 41, 21
    grid = [np.linspace(data[:, i].min(), data[:, i].max(), g + 1)
            for i, g in enumerate(Ngrid)]

    H, (xbins, ybins) = np.histogramdd(data, bins=grid)

    # The bin order in the histogram now becomes z, y, x!!!
    H = H.T

    # Now we try to recover, and predict
    expected_mass = np.zeros((test.size, ),
                             dtype=[('MASS', '>f4'), ('MASS_err', '>f4',
                                                      (2, ))])

    for j, s in enumerate(test['LOSVD_dist']):
        # using the distribution of sigmas we make P(s) and ds
        Ps, ds = np.histogram(np.log10(s), bins=50, density=True)
        # and get new new sigmas from the histogram2d
        centers = (ds[1:] + ds[:-1]) / 2.0

        # now we have to make P(mt|s) -- from the training sample
        # will need to do this a bunch of times for each sigma in centers
        iy = np.digitize(centers, ybins)  #change this line
        # now we make P(mt) using P(mt|s) and P(s)ds from above
        Pmt = np.zeros(xbins.size - 1)
        Psds = Ps * np.diff(ds)
        for i in range(iy.size):
            try:
                Pmt_s = H[iy[i]] / H[iy[i]].sum()
            except IndexError:
                Pmt_s = np.zeros(xbins.size - 1)
            Pmt += Pmt_s * Psds[i]

        # now we can calculate the expected mass
        centers = (xbins[:-1] + xbins[1:]) / 2.
        norm = np.sum(Pmt * np.diff(xbins))
        M_expect = np.sum(centers * Pmt * np.diff(xbins)) / norm

        variance = np.sum(
            (centers - M_expect)**2 * Pmt * np.diff(xbins)) / norm

        expected_mass['MASS'][j] = M_expect
        expected_mass['MASS_err'][j] = [M_expect - np.sqrt(variance),
                                        M_expect + np.sqrt(variance)]

        #mask = np.where(np.isnan(expected_mass['MASS']))[0]
        #expected_mass = np.delete(expected_mass, mask)
        #test = np.delete(test, mask)

    return expected_mass


if __name__ == "__main__":
    prob1d()
