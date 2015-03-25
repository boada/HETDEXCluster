import glob
import pandas as pd
import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca

c = 2.99E5 # speed of light in km/s

def findSeperationSpatial(data, center):
    ''' Finds the distance to all of the galaxies from the center of the
    cluster in the spatial plane. Returns values in Mpc.

    '''

    # Add a new column to the dataframe
    data['seperation'] = 0.0
    for row in data.iterrows():
        sepDeg = aco.calcAngSepDeg(center[0], center[1], row[1]['ra'],
                row[1]['dec'])
        sepMpc = sepDeg * aca.da(row[1]['redshift'])/57.2957795131
        data['seperation'][row[0]] = sepMpc

    return data

def findClusterCenterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator.

    '''
    x = np.copy(data['redshift'].values)
    return ast.biweightLocation(x, tuningConstant=6.0)


def findLOSV(data):
    ''' Finds the line of sight velocity for each of the galaxies.

    '''

    c = 2.99E5 # speed of light in km/s

    avgz = findClusterCenterRedshift(data)

    # Add a new column to the dataframe
    data['LOSV'] = 0.0
    for row in data.iterrows():
        data['LOSV'][row[0]] = c *(row[1]['redshift'] - avgz)/(1 + avgz)

    return data

def split_list(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list will have an extra element.

    '''
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
        for i in range(wanted_parts) ]

def rejectInterlopers(data):
    ''' Does all of the work to figure out which galaxies don't belong. Makes
    several sorted copies of the dataframe and then applies the fixed gapper
    method.

    '''

    # make some copies so we can sort them around
    sepSorted = data.sort('seperation', ascending=True)
    # How many parts to break into
    parts = len(data)//15
    splitData = split_list(data, parts)

    # Now we sort the parts by LOSV and find the rejects
    interlopers = []
    for part in splitData:
        # sort by LOSV
        LOSVsorted = part.sort('LOSV', ascending=True)
        rejected = True
        while rejected:
            # Find the difference between all of the neighboring elements
            difference = np.diff(LOSVsorted['LOSV'])
            # If diff > 1000 reject
            rejects = abs(difference) > 1000
            # Now remove those items
            indices = np.where(rejects == True)
            #print LOSVsorted['LOSV']
            #print difference
            #print indices[0]
            if rejects.any() == True:
                # Always take the more extreme index
                for index, i in enumerate(indices[0]):
                    if (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                        abs(LOSVsorted['LOSV'][LOSVsorted.index[i+1]])) > 0:
                            pass
                    elif (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                        abs(LOSVsorted['LOSV'][LOSVsorted.index[i+1]])) < 0:
                            indices[0][index] = i+1

                #print LOSVsorted.index[list(indices[0])]
                dataframeIndex = list(LOSVsorted.index[list(indices[0])])
                LOSVsorted = LOSVsorted.drop(dataframeIndex)
                interlopers += dataframeIndex
            else:
                rejected = False
    print 'interlopers',interlopers
    return data.drop(interlopers)

def rejectInterlopers_group(data, sigmav=500)
    ''' This method is given in Wilman 2005 and Connelly 2012 and talked about
    in the draft of the cluster paper. It doesn't look like it is full fleshed
    out.

    '''

    deltaZmax = 2 * simgav / c
    avgz = findClusterCenterRedshift(data)

    deltaRmax = (c * deltaZmax)/(10*(1 + avgz)*aca.H0*aca.Ez(avgz)) # 1/Mpc

    deltaThetamax = 206265 * deltaRmax * aca.da(avgz) # arcseconds





center = 191.1050125, 16.7966666667

seperated = findSeperationSpatial(matched, center)
losv = findLOSV(seperated)

cleaned = rejectInterlopers(losv)

