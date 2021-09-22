from Preliminary1D import InputGeometry
from Preliminary1D import InputVariables
from Preliminary1D import InputSoluteParameters
import numpy as np
from pytest import approx
from Solver1D import rad_solver_1D



def test_Solver1D():

    """
    # Test conservation of r-a-d solver (no production or degradation of solute)
    
    """

    #--  setup problem

    geom_dict = InputGeometry(L = 1e-3, T = 1e3, nx = 1000)

    n, w, uw, parameters_dict = InputVariables(geom_dict)

    c0, param_dict = InputSoluteParameters(parameters_dict,  alpha = 0, kappa = 0, delta = 0)

    #-- initial conditions, variables 
    nx = param_dict['nx']
    c0 = 0.1*np.random.rand(nx)

    #-- calculate total initial solute
    cint = np.sum(c0)

    #-- solve numerically using matrix solver
    nsteps = 10

    for t in range(nsteps):
        c = rad_solver_1D(c0, n, w, uw, param_dict)
        c0 = c.copy()

    #-- calculate total final solute

    cfin = np.sum(c) 

    #-- check conservation

    assert cfin == approx(cint, rel = 1e-3)

    return