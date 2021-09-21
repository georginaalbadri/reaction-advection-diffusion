import numpy as np
from pytest import approx
from Preliminary1D import InputSoluteParameters
from Preliminary1D import InputVariables



def test_InputVariables():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    assert len(parameter_dict) == 7
    assert list(parameter_dict.keys()) == ['L', 'T', 'nx', 'dx', 'dx2', 'phi', 'm']
    assert np.shape(n) == (nx, )
    assert np.shape(w) == (nx, )
    assert np.shape(uw) == (nx, )

    return




def test_ImportSoluteParameters():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

    assert len(param_dict) == 14
    assert list(param_dict.keys()) == ['L', 'T', 'nx', 'dx', 'dx2', 'phi', 'm', 'dt', 'c0', 'D', 'alpha', 'kappa', 'K', 'delta']
    assert np.shape(c0) == (nx, )

    return