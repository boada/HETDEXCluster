import numpy as np
from scipy.optimize import curve_fit
from astLib.astStats import bootstrap
from astLib import astCalc as aca

def correctFit(x, a):
    ''' This is the general form of the VD-mass relation that we are going to
    use. This function is called by the curve fitter to correct the A1D
    parameter for our specific data case.

    '''

    return a * (x/1e15)**(1/3.)

def wrapper_curveFit(y, x, a):
    ''' This is a simple wrapper which lets us figure out the bootstrap CI for
    the parameter given from correctFit.

    '''

    opt_parms, parm_cov = curve_fit(correctFit, x, y, p0=[a])
    return opt_parms[0]

def mass(vd, z, a):

    m200 = (vd/a)**(3) * 1e15/aca.Ez(z)
    return m200
