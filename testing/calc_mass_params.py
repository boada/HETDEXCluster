import numpy as np
from scipy.optimize import curve_fit
from astLib.astStats import bootstrap
from astLib import astCalc as aca
import h5py as hdf


# aardvark simulation cosmology
aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def mass(s1d, clusz, a1d=1082, a=1/3.):
    ''' This is the general form of the VD-mass relation that we are going to
    use. This function is called by the curve fitter to correct the A1D
    parameter for our specific data case.

    '''

    # THIS IS NOT TO BE USED for THE ESTIMATION!

    hz = aca.H0 * aca.Ez(clusz)/100.

    return 1e15/hz * (s1d/a1d) ** 1./a

def vd(data, a1d=1082, a=1/3.):
    ''' This is the general form of the VD-mass relation that we are going to
    use. This function is called by the curve fitter to correct the A1D
    parameter for our specific data case.

    '''

    #mass = data[:,0]
    #clusz = data[:,1]
    mass = data['M200c']
    clusz = data['ZSPEC']

    hz = np.array([aca.H0 * aca.Ez(z)/100 for z in clusz])

    return  a1d * ((hz * mass)/1e15)**a


def wrapper_curveFit(mass, s1d):
    ''' This is a simple wrapper which lets us figure out the bootstrap CI for
    the parameter given from correctFit.

    '''

    opt_parms, parm_cov = curve_fit(mass, mass, s1d, p0=[1082, 1/3.] args=(clusz))
    return opt_parms[0]

