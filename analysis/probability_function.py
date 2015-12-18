import numpy as np

class probabilityFactory:
    def __init__(self, *args):
        ''' The first column should correspond to the true mass. After that,
        the order is y, z, ...

        '''

        print 'I will use', len(args) - 1, 'training features'

        self.numFeatures = len(args) - 1
        self.trainingData = np.column_stack(args)

    def train(self, Ngrid):
        ''' The first column has to correspond to the true mass. After that,
        the order is y, z, ...

        '''

        if not len(Ngrid) == numFeatures:
            raise Exception('Number of Features does not match grid')

        grid = [np.linspace(data[:,i].min(), data[:,i].max(), g+1) for i,g in
                    enumerate(Ngrid)]
        self.H, self.bins = np.histogramdd(self.trainingData, bins=grid)
        # The bin order in the histogram now becomes z, y, x!!!
        self.H = self.H.T

    def predict(self, *args):
        if not len(args) == numFeatures:
            raise Exception('Number of Features does not match grid')
        # using the distribution of sigmas we make P(s) and ds
        Ps, ds = np.histogram(np.log10(args[1]), bins=50, density=True)
        centers = (ds[1:] + ds[:-1])/2.0

        # now we have to make P(m|s,z) -- from the training sample
        # will need to do this a bunch of times for each sigma in centers
        for i, f in enumerate(zip(args, bins), start=1):
            iy  = np.digitize(centers, self.ybins)
            iz = np.digitize([z], self.zbins)[0]
        



        # now we make P(m) using P(m|s,z) and P(s)ds from above
        Pm = np.zeros(self.xbins.size -1)
        Psds = Ps*np.diff(ds)
        for i, b in enumerate(iy):
                try:
                    if self.H[iz-1, b-1, :].sum() == 0.0:
                        Pm_sz = np.zeros(self.xbins.size -1)
                    else:
                        Pm_sz = self.H[iz-1, b-1, :]/ self.H[iz-1, b-1, :].sum()
                except IndexError:
                    Pm_sz = np.zeros(self.xbins.size -1)
                #print Pm_sz
                Pm += Pm_sz * Psds[i]

        # now we can calculate the expected mass
        centers = (self.xbins[:-1] + self.xbins[1:])/2.
        norm = np.sum(Pm * np.diff(self.xbins))
        M_expect = np.sum(centers * Pm * np.diff(self.xbins))/norm


if __name__ == "__main__":

