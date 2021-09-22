import numpy as np
from pytest import approx
from Preliminary1D import InputGeometry
from Preliminary1D import Nondimensionalise
from Preliminary1D import InputSoluteParameters
from Preliminary1D import InputVariables


def test_InputGeometry():

    #test case
    geom_dict = InputGeometry()

    assert len(geom_dict) == 5
    assert list(geom_dict.keys()) == ['L', 'T', 'nx', 'dx', 'dx2']

    return 


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



def test_Nondimensionalise():

    #-- set problem geometry 
    L = 5e-3  #length-scale [m]
    T = 1e3  #timescale [s]
    nx = 500 #number of grid points
    dx = 1/nx #dimensionless grid spacing
    dx2 = dx*dx
    geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    #test cases
    dimval = 1
    cM = 1e-9

    #function results
    val1 = Nondimensionalise(parameter_dict, cM = cM, dimval = dimval, param_name = 'D')
    val2 = Nondimensionalise(parameter_dict, cM = cM, dimval = dimval, param_name = 'alpha')
    val3 = Nondimensionalise(parameter_dict, cM = cM, dimval = dimval, param_name = 'delta')
    val4 = Nondimensionalise(parameter_dict, cM = cM, dimval = dimval, param_name = 'K')
    
    #manual results 
    D = dimval * T / L**2
    alpha = dimval * T / cM
    delta = dimval * T
    K = dimval / cM

    #compare function and manual results
    assert val1 == approx(D)
    assert val2 == approx(alpha)
    assert val3 == approx(delta)
    assert val4 == approx(K)

    return 