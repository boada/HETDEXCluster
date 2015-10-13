import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
from numpy.lib import recfunctions as rfns
from sklearn import mixture

# aardvark simulation cosmology
aca.H0 = 72
aca.OMEGA_M0 = 0.23
aca.OMEGA_L0 = 0.77

def findLOSVD(data):
    if data.size >= 15:
        data['LOSVD'] = ast.biweightScale_test(data['LOSV'],
                tuningConstant=9.0)
        try:
            data['LOSVD_err'] = ast.bootstrap(data['LOSV'],
                    ast.biweightScale_test, resamples=10000, alpha=0.32,
                    output='errorbar', tuningConstant=9.0)
        except ZeroDivisionError:
            data['LOSVD_err'] = [0.0, 0.0]
    elif data.size >=5:
        data['LOSVD'] = ast.gapperEstimator(data['LOSV'])
        data['LOSVD_err'] = ast.bootstrap(data['LOSV'], ast.gapperEstimator,
                resamples=10000, alpha=0.32, output='errorbar')
    else:
        data['LOSVD'] = 0.0
        data['LOSVD_err'] = [0.0, 0.0]

    return data

def findLOSVDgmm(data):

    LOSV = data['LOSV']
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 4)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    #cv_types = ['diag']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a mixture of Gaussians with EM
            gmm = mixture.GMM(n_components=n_components,
                    covariance_type=cv_type)
            gmm.fit(LOSV)
            bic.append(gmm.bic(LOSV))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm

    # figure things out -- this comes from the wikipedia article on mixture
    # distributions. I have no idea if it is completely trustworthy.
    #covars = best_gmm.covars_.ravel()
    #weights = best_gmm.weights_.ravel()
    #means = best_gmm.means_.ravel()
    #wmeans = np.sum(weights*means)

    #parts = weights * ((means - wmeans)**2 + covars)
    #newLOSVD = np.sqrt(np.sum(parts))

    #data['R200'] = best_gmm.n_components

    ## now we resample and then see
    dx = np.linspace(LOSV.min()-100,LOSV.max()+100,1000)
    logprob, responsibilities = best_gmm.score_samples(dx)
    pdf = np.exp(logprob)

    normedPDF = pdf/np.sum(pdf)

    u = np.sum(dx*normedPDF)
    data['LOSVDgmm'] = np.sqrt(np.sum(normedPDF*(dx-u)**2))
    #data['LOSVD'] = newLOSVD
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
    data['R200'] = np.sqrt(3) * (LOSVD)/(10*aca.H0 * aca.Ez(avgz))
    return data

def splitList(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list wil have an extra element.

    '''

    length = len(alist)
    return [alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i
        in range(wanted_parts)]


def z2v(z, zc):
    """Convert the redshift to km/s relative to the cluster center"""
    return 2.99792458e5*(z-zc)/(1+zc)


def shifty_gapper(r, z, zc, vlimit=10000, ngap=30, glimit=1000):
   """Determine cluster membersip according to the shifty
      gapper technique
      The shifty gapper technique includes the following steps:
      1.  Remove outliers outside of +- vlimit
      2.  Locate the Ngap targets with radius close to r_i
      3.  Within this sample of N targets, identify sources with
          |v_pec| < |v_pec_i|
      4.  Meaure the largest gap within this subset of targets
      5.  If v_gap is larger than glimit, reject the source
      Parameters
      -----------
      vlimit: float
         Initial limit in velocity offset from the cluster center
      ngap: int
         Number of sources to use in in estimating the gap
      glimit:  float
         Maximum gap size for rejecting objects

      Parameters
      -----------
      incluster: ndarray
          Returns a boolean array where objects in the cluster have
          a value of True
   """

   #convert to the velocity scale
   v = z2v(z,zc)

   #limit the source to only sources within the vlimit
   vmask = abs(v) < vlimit

   nobj=len(r)
   incluster=np.zeros(nobj, dtype=bool)

   if nobj<ngap:
      raise Exception('Number of sources is less thant number of gap sources')

   for i in range(nobj):
     if abs(v[i])<vlimit:
       #find the ngap closest sources
       r_j=abs(r[vmask]-r[i]).argsort()
       vg=v[vmask][r_j[0:ngap]]

       #find the sources with |v_pec| < |v_pec_i|
       mask=abs(vg)<=abs(v[i])
       if mask.sum()>1:
          vg=vg[mask]
          #now sort these sources and find the biggest gap
          vg.sort()
          if np.diff(vg).max()<glimit: incluster[i]=True
       else:
          incluster[i]=True

   return incluster


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
    data = rfns.append_fields(data, ['SEP', 'CLUSZ', 'LOSV', 'LOSVD',
        'LOSVDgmm', 'MASS', 'R200', 'NGAL'], [newData, newData, newData,
            newData, newData, newData, newData, newData], dtypes='>f4',
        usemask=False)

    newnewData = np.zeros(data.size, dtype=[('LOSVD_err', '>f4', (2,)),
        ('LOSVDgmm_err', '>f4', (2,))])
    data = rfns.merge_arrays((data, newnewData), usemask=False,
            asrecarray=False, flatten=True)

    return data
