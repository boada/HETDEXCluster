import numpy as np
from astLib import astCoords as aco
def easy_nn(X, n=7):
    N = X.shape[0]
    neighbors = np.zeros((N,n), dtype=int)
    distances = np.zeros((N,n), dtype=float)
    for i in range(N):
        #d = aco.calcAngSepDeg(X['RA'][i], X['DEC'][i], X['RA'], X['DEC'])
        d = np.sqrt(np.sum((X[i] - X) **2, axis=1 ))
        s = np.argsort(d)
        neighbors[i] = s[:n]
        distances[i] = d[s[:n]]
    return neighbors, distances

def vectorized_nn(X, n=7):
    distances = np.zeros(X.shape[0])
    XXT = np.dot(X, X.T)
    Xii = XXT.diagonal()
    D = Xii - 2 * XXT + Xii[:, np.newaxis]
    s = np.argsort(D, axis=1)[:,n]
    distances = D[s]
    return distances
    

