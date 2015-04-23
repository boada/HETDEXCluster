import numpy as np
import h5py as hdf
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
from numpy.lib import recfunctions as rfns

# aardvark simulation cosmology
aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def calcVD_big(data):
    ''' Calculates the velocity dispersion for a large group of members '''

    return ast.biweightScale(data, tuningConstant=9.0)

def calcVD_small(data):
    ''' Calculates the velocity dispersion for a small group of members '''

    return ast.gapperEstimator(data)

def calcVD_test(data):
    return ast.biweightScale_test(data, tuningConstant=9.0)

def calc_mass_Saro(data):
    ''' Calculates the mass using the scaling relation from Saro2013. '''

    avgz = data['CLUSZ'][0]
    vd = data['LOSVD'][0]

    a = 939
    b = 2.91
    c = 0.33

    m200 = (vd**b/(a * ((aca.H0 * aca.Ez(avgz))/72)**c)**b)
    return m200 * 1e15

def calc_mass(data):
    ''' Calculates the dynamical mass (m200) and radius (r200) of the galaxy
    cluster. I suppose it also calculates the velocity dispersion.

    '''

    if len(data) > 10:
        vd = calcVD_big(data['LOSV'])

    else:
        vd = calcVD_small(data['LOSV'])

    avgz = data['CLUSZ'][0]

    r200 = np.sqrt(3) * vd /(10*aca.H0 * aca.Ez(avgz))
    a = 3 * np.sqrt(3) * 1000**3 * 3.08E19/(10*aca.H0 * aca.Ez(avgz) * \
            6.67384E-11)

    m200 = a * vd**3

    return vd, m200/1.9891E30, r200

def findSeperationSpatial(data, center, unit='Mpc'):
    ''' Finds the distance to all of the galaxies from the center of the
    cluster in the spatial plane. Returns the values in Mpc.

    '''

    data['SEP'] = aco.calcAngSepDeg(center[0], center[1], data['RA'],
            data['DEC'])

    if unit == 'Mpc':
        for index, value in enumerate(data['Z']):
            data['SEP'][index] *= aca.da(value)/57.2957795131
    elif unit == 'arcsecond':
        data['SEP'] *= 3600

    return data

def findClusterCenterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator. Puts the result in the ClusZ column of the data
    array.

    '''

    data['CLUSZ'] = ast.biweightLocation(data['Z'], tuningConstant=6.0)
    #return ast.biweightClipped(data['Z'], 6.0, 3)['biweightLocation']
    #return ast.biweightLocation(data['Z'], tuningConstant=6.0)

    return data

def findLOSV(data, ClusZ=None):
    ''' Finds the line of sight velocity for each of the galaxies and puts it
    in the LOSV column of the data array.

    '''

    c = 2.99e5 # speed of light in km/s
    if ClusZ == None:
        data['LOSV'] = c * (data['Z'] - data['CLUSZ'])/(1 + data['CLUSZ'])
        return data
    else:
        losv = c * (data['Z'] - ClusZ)/(1 + ClusZ)
        return losv

def splitList(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list wil have an extra element.

    '''

    length = len(alist)
    return [alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i
        in range(wanted_parts)]

def rejectInterlopers(data):
    ''' Does all of the work to figure out which galaxies don't belong. Makes a
    couple sorted copies of the data and then applies the fixed gapper method.

    '''

    # sort the part that we want in place in the array
    data = np.sort(data, order='SEP')
    parts = len(data)//15
    splitData = splitList(data, parts)

    for part in splitData:
        # sort by LOSV
        LOSVsorted = np.sort(part, order='LOSV')
        rejected = True
        while rejected:
            # Find the difference between all the neighboring elements
            difference = np.diff(LOSVsorted['LOSV'])
            # If diff > 1000 km/s, reject
            rejects = abs(difference) > 1000
            # Now remove those items
            indices = np.where(rejects)
            if rejects.any():
                # Always take the more extreme index
                for index, i in enumerate(indices[0]):
                    if (abs(LOSVsorted['LOSV'][i]) -
                            abs(LOSVsorted['LOSV'][i+1]) > 0):
                        pass
                    elif (abs(LOSVsorted['LOSV'][i]) -
                            abs(LOSVsorted['LOSV'][i+1]) < 0):
                        indices[0][index] = i + 1

                print indices[0]
                LOSVsorted = np.delete(LOSVsorted, indices[0])
            else:
                rejected = False
                try:
                    result = np.append(result, LOSVsorted)
                except NameError:
                    result = LOSVsorted
    return result

def rejectInterlopers_group(data, sigmav=500):
    ''' This method is given in Wilman+2005 and Connelly+2012 and talked about
    in the draft o fhte cluster paper. It should reject interlopers when there
    are fewer than 15 members and we can't use the other method.

    '''

    c = 2.99e5

    avgz = data['CLUSZ'][0]

    deltaZmax = 2 * sigmav/ c
    deltaRmax = (c *deltaZmax)/(3.5*(1+avgz)*aca.H0*aca.Ez(avgz)) # 1/Mpc
    deltaThetaMax = 206265 * deltaRmax * aca.da(avgz) # arcseconds

    # Apply limits
    selected = (abs(data['Z'] - avgz) <= deltaZmax) & (data['SEP'] <=\
            deltaThetaMax)

    return data[selected]

def updateArray(data):
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    l = len(data)
    data = rfns.append_fields(data, ['SEP', 'CLUSZ', 'LOSV'], [np.ones(l) * -1,
        np.ones(l) * -1, np.ones(l) * -1], usemask=False)
    return data

def updateArray2(data):
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    newData = -np.ones(len(data))
    data = rfns.append_fields(data, ['CRA', 'CDEC', 'CZ', 'VRMS',
        'NGALS', 'M200', 'R200', 'CLUSZ', 'LOSV', 'LOSVD', 'MASS'], [newData,
            newData, newData, newData, newData, newData, newData, newData,
            newData, newData, newData], dtypes='>f4', usemask=False)

    return data

def countGalaxies(data):
    ''' Counts the number of galaxies we have associated with each of the
    halos.

    '''
    from collections import Counter
    count = Counter(data['HALOID'])
    return count


