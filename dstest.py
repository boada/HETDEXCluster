import numpy

def DSmetric():
    """ Calculates the delta squared values as given in equation 1 of Dressler
    et al. 1988. This uses both the position and velocities to give a measure
    of substructure by identifying galaxies with do not follow the cluster
    velocity distribution.

    @type x: list
    @param x: x (RA) position, accepts a list of RA values in decimal degrees
    @type y: list
    @param y: y (DEC) position, accepts a list of DEC values in decimal degrees
    @type v: list
    @param v: velocities associated with the galaxies
    @rtype: list
    @return: a list of delta squared values

    """




def DStest():
    """ Computes the delta deviation for an entire cluster. See Dressler et al.
    1988 for a complete description. A result close to the number of galaxies
    present in the cluster indicates no substructure.

    """

    dsm = DSmetric()

