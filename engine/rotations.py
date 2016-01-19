import numpy as np
import random

def randomR():
    ''' randomly choose a location on the unit hemisphere, r = [x,y,z] with
    -1<x<1, -1<y<1, and 0<z<1 and xˆ2+yˆ2+zˆ2 = 1.

    '''

    phi = 2*np.pi*random.random()
    theta = np.arccos(np.random.random())
    r = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
    return r

def normVector(v):
    return v/ np.sqrt(np.dot(v,v))


def planeVectors(v):
    ''' returns two 3d vectors that are orthogonal to the vector v, orthogonal
    to each other, and normalized.

    '''

    if v[0] == 1. or v[1] == 1. or v[2] == 1.:
        # v lies along a box direction
        r1 = normVector([v[1], v[2], v[0]])
    else:
        r1 = normVector([v[0], -v[1], 0])
    r2 = normVector(np.cross(v, r1))
    return [r1, r2]
