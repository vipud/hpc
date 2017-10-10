import numpy as np
from timeit import default_timer as timer
from Ihelp import *
from numba import cuda
from math import exp, sqrt


def driver(MCUDA, pinned=False):
    paths = np.zeros((Number_of_Paths, Number_of_Steps + 1), order='F')
    paths[:, 0] = SP
    Delta_T = Maturity / Number_of_Steps

    if pinned:
        with cuda.pinned(paths):
            time_s = timer()
            MCUDA(paths, Delta_T, IR, Beta)
            time_a = timer()
    else:
        time_s = timer()
        MCUDA(paths, Delta_T, IR, Beta)
        time_a = timer()

    stk = paths[:, -1]
    pOff = np.maximum(paths[:, -1] - K, 0)
    o_Price = np.mean(pOff)*exp(Maturity * -IR)
    print('error ', np.std(stk)/ sqrt(Number_of_Paths))
    print('payoff ', np.mean(pOff))
    print('Option Pice', o_Price)
    print('Run Time')
    print(time_a - time_s)




