import numpy as np
from calc_cluster_props import *

number = 1000

result = np.zeros((number-5), dtype=[('len', int), ('vduc', float), ('vdc',
    float), ('vdg', float), ('vdbi',float), ('vdgap', float), ('vdt', float)])

#vd**2 = s**2 pi**2 /3 logistic
#vd**2 = s**2 pi**2/6 gumbel
#vd = sqrt(lam) poisson

for idx, i in enumerate(range(5,number)):

    #losv = np.random.normal(loc=2, scale=0.5, size=i)
    #losv = np.random.logistic(loc=2, scale=0.5, size=i)
    #losv = np.random.gumbel(loc=2, scale=0.5, size=i)
    #losv = np.random.poisson(lam=9, size=i)
    losv = np.random.laplace(loc=2, scale=0.5, size=i)

    try:
        result[idx]['len'] = i
        result[idx]['vduc'] = np.std(losv)
        result[idx]['vdc'] = np.std(losv, ddof=1)
        result[idx]['vdg'] = np.std(losv, ddof=1.5)
        result[idx]['vdbi'] = calcVD_big(losv)
        result[idx]['vdgap'] = calcVD_small(losv)
        result[idx]['vdt'] = calcVD_test(losv)
    except:
        pass
