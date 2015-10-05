import numpy as np
from astLib import astCoords as aco
def easy_nn(X, n=7):
    N, D = X.size, X.size
    neighbors = np.zeros((N,n), dtype=int)
    distances = np.zeros((N,n), dtype=float)
    for i in range(N):
        d = aco.calcAngSepDeg(X['RA'][i], X['DEC'][i], X['RA'], X['DEC'])
        s = np.argsort(d)
        neighbors[i] = s[:n]
        distances[i] = d[s[:n]]
    return neighbors, distances

