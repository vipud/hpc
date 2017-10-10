import numpy as np
from math import exp, sqrt
from numba import cuda, vectorize
from pyculib import rand as curand
from Driver import driver


class Run(object):
    def __init__(self):
        self.path_s = driver(MCUDA)

    def ree(self):
        return self.path_s


@vectorize(['f8(f8, f8, f8, f8, f8)'], target='cuda')
def stepper(last, delta, c0, c1, alpha):
    return last * exp(c0*delta+c1*alpha)


def MCUDA(paths, delta, IR, Beta):
    n = paths.shape[0]
    prng = curand.PRNG(rndtype=curand.PRNG.MRG32K3A)
    d_norm = cuda.device_array(n, dtype=np.double)
    c0 = IR - .5 * Beta **2
    c1 = Beta * sqrt(delta)
    d_previous = cuda.to_device(paths[:, 0])
    for i in range(1, paths.shape[1]):
        prng.normal(d_norm, mean=0, sigma=1)
        d_path = cuda.to_device(paths[:, i])
        stepper(d_previous, delta, c0, c1, d_norm, out=d_path)
        d_path.copy_to_host((paths[:,i]))
        d_previous = d_path



