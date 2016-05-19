import numpy as np

def convert(m, merr, b, berr, lam):
    llam = np.log10(lam)

    step1 = m * llam
    step1_err = merr * llam

    step2 = b + step1
    step2_err = np.sqrt(berr**2 + step1_err**2)


    step3 = 10**step2
    step3_err = 2.303 * step3 * step2_err

    return step3/1e14, step3_err/1e14
