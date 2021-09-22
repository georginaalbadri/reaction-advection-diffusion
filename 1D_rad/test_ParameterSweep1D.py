import numpy as np
from pytest import approx
from Preliminary1D import InputGeometry, InputVariables, InputSoluteParameters
from ParameterSweep1D import ParamSweep1D


def test_ParamSweep1D():

    geometry_dict = InputGeometry(L = 5e-3, T = 1e3, nx = 200)

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

    #-- test case
    param_list = [1e-13, 1e-12, 1e-11]
    Tmax = 0.1
    dt = param_dict['dt']
    nsteps = int(Tmax / dt)


    clist = ParamSweep1D(c0, n, w, uw, param_list, param_dict, param_name = 'alpha', nsteps = nsteps)

    assert len(clist) == len(param_list)
    assert len(clist[0]) == param_dict['nx']

    return
