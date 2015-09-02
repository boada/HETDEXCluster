import multiprocessing
import numpy as np
import h5py as hdf
from os import environ

def mp_worker():
    pass

def mp_worker_wrapper(args):
    return mp_worker(*args)

if __name__ == "__main__":
    truthPath = './../data/buzzard_v1.0/allbands/truth/'
    haloPath = './../data/buzzard_v1.0/allbands/halos/'

    p = multiprocessing.Pool(multiprocessing.cpu_count())
