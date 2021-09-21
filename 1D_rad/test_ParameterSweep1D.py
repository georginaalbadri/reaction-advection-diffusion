import numpy as np
from pytest import approx
from Preliminary1D import InputVariables
from Preliminary1D import InputSoluteParameters
from ParameterSweep1D import ParamSweep1D, ParamSweepPlot1D


def test_ParamSweep1D():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

    #-- test case
    param_list = [1e-13, 1e-12, 1e-11]
    Tmax = 0.1
    dt = param_dict['dt']
    nsteps = int(Tmax / dt)


    clist = ParamSweep1D(c0, n, w, uw, param_list, param_dict, param_name = 'alpha', nsteps = nsteps)

    assert len(clist) == len(param_list)
    assert len(clist[0]) == nx

    return
