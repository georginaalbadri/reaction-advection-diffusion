import numpy as np
from SensitivityAnalysis1D import multiprocessing_func


def test_multiprocessing_func():

    # test  case
    D = 1e-10
    alpha = 1e-11
    kappa = 1e-14
    K = 0.1
    delta = 1e-5 

    param_values = [D, alpha, kappa, K, delta]

    ctot = multiprocessing_func(param_values)

    assert np.isnan(ctot) == False 

    return