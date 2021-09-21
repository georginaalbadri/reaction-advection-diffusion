import numpy as np
import matplotlib.pyplot as plt
from pytest import approx
from Solver1D import rad_solver_1D
import sys



def test_Solver1D():

    """
    # Test conservation of r-a-d solver (no production or degradation of solute)
    
    """

    #--  parameter values
    ny = 1000
    dy = 1/ny
    dy2 = dy*dy
    dt = 0.001

    alpha = 0
    delta = 0
    kappa = 0
    K = 0.1

    D = 0.1

    #-- initial conditions, variables 
    c0 = 0.1*np.random.rand(ny)
    m = 0.03
    phi = 1 - m
    n = 0.03 * np.ones((ny))
    w = phi - n 
    vw = np.zeros((ny))

    param_dict = {'L':1, 'T':1e3, 'ny':ny, 'dy':dy, 'dy2':dy2, 'phi':phi, 'm':m, 'dt':dt, \
        'cint':c0, 'D':D, 'alpha':alpha, 'kappa':kappa, 'K':K, 'delta':delta}

    #-- calculate total initial solute
    cint = np.sum(c0)

    #-- solve numerically using matrix solver
    nsteps = 10

    for t in range(nsteps):
        c = rad_solver_1D(c0, n, w, vw, param_dict)
        c0 = c.copy()

    #-- calculate total final solute

    cfin = np.sum(c) 

    #-- check conservation

    assert cfin == approx(cint, rel = 1e-3)

    return