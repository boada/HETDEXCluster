import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
from numpy.lib import recfunctions as rfns

# aardvark simulation cosmology
aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def findLOSVD(data):
    if data.size >= 15:
        data['LOSVD'] = ast.biweightScale_test(data['LOSV'],
                tuningConstant=9.0)
    elif data.size >=5:
        data['LOSVD'] = ast.gapperEstimator(data['LOSV'])
    else:
        data['LOSVD'] = 0.0

    return data

def calc_mass_Saro(data):
    ''' Calculates the mass using the scaling relation from Saro2013. '''

    avgz = data['CLUSZ'][0]
    vd = data['LOSVD'][0]

    a = 939
    b = 2.91
    c = 0.33

    m200 = (vd**b/(a * ((aca.H0 * aca.Ez(avgz))/72)**c)**b)
    return m200 * 1e15

def calc_mass_Evrard(data, A1D = 1082.9, alpha=0.3361):
    ''' Calculates the mass using the scaling relation from Saro2013. Really it
    is using the relationship from Munari2013.'''

    avgz = data['CLUSZ'][0]
    vd = data['LOSVD'][0]
    if vd == 0.0:
        data['MASS'] = 0.0
        return data
    data['MASS'] = 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)
    return data

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

def findClusterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator. Puts the result in the ClusZ column of the data
    array.

    '''
    if data.size >= 5:
        data['CLUSZ'] = ast.biweightLocation(data['Z'], tuningConstant=6.0)
    else:
        data['CLUSZ'] = 0.0
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

def findR200(data):
    LOSVD = data['LOSVD'][0]
    avgz = data['CLUSZ'][0]
    return np.sqrt(3) * (LOSVD)/(10*aca.H0 * aca.Ez(avgz))

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

                #print indices[0]
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

    deltaZmax = 2 * sigmav * (1 + avgz) / c
    #deltaRmax = (c *deltaZmax)/(3.5*(1+avgz)*aca.H0*aca.Ez(avgz)) # 1/Mpc
    deltaRmax = (c * deltaZmax)/(9.5*aca.H0*aca.Ez(avgz)) # Mpc

    # Apply limits
    selected = (abs(data['Z'] - avgz) <= deltaZmax) & (data['SEP'] <=\
            deltaRmax)

    return data[selected]

def updateArray(data):
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    newData = -np.ones(data.size)
    data = rfns.append_fields(data, ['SEP', 'CLUSZ', 'LOSV', 'LOSVD', 'MASS'],
            [newData, newData, newData, newData, newData], dtypes='>f4',
            usemask=False)
    return data
def updateArray2(data):
    ''' Makes the new fields that we are going to add things into in the
    functions above. This should only be called once.

    '''

    newData = -np.ones(data.size)
    data = rfns.append_fields(data, ['IDX', 'SEP', 'CLUSZ'],
            [newData, newData, newData], dtypes='>f4',
            usemask=False)
    return data 
