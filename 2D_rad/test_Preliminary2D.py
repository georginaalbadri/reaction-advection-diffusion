import numpy as np
from pytest import approx
from Preliminary2D import InputSoluteParameters
from Preliminary2D import InputVariables



def test_InputVariables():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    ny = nx
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'ny':ny, 'dx':dx, 'dx2':dx2} 

    n, w, uw, vw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    assert len(parameter_dict) == 8
    assert list(parameter_dict.keys()) == ['L', 'T', 'nx', 'ny', 'dx', 'dx2', 'phi', 'm']
    assert np.shape(n) == (nx, ny)
    assert np.shape(w) == (nx, ny)
    assert np.shape(uw) == (nx, ny)
    assert np.shape(vw) == (nx, ny)

    return




def test_ImportSoluteParameters():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    ny = nx
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'ny':ny, 'dx':dx, 'dx2':dx2} 

    n, w, uw, vw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

    assert len(param_dict) == 15
    assert list(param_dict.keys()) == ['L', 'T', 'nx', 'ny', 'dx', 'dx2', 'phi', 'm', 'dt', 'c0', 'D', 'alpha', 'kappa', 'K', 'delta']
    assert np.shape(c0) == (nx, ny)

    return