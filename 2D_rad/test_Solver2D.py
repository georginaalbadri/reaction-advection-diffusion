import numpy as np
import matplotlib.pyplot as plt
from pytest import approx
from Solver2D import rad_solver_2Daxi




def test_Solver2D():

    """
    # Test conservation of 2D axisymmetric r-a-d solver (no production or degradation of solute)
    
    """

    #--  parameter values
    ny = 100
    dy = 0.1
    dy2 = dy*dy
    dt = 0.1
    nx = ny 

    alpha = 0
    delta = 0
    kappa = 0
    K = 0.1

    D = 0.1

    #-- initial conditions, variables 
    c0 = 0.1*np.random.rand(ny, ny)
    m = 0.03
    phi = 1 - m
    n = 0.03 * np.ones((ny, ny))
    w = phi - n 
    uw = np.zeros((ny, ny))
    vw = np.zeros((ny, ny))

    param_dict = {'L':1, 'T':1e3, 'nx':nx, 'ny':ny, 'dy':dy, 'dy2':dy2, 'phi':phi, 'm':m, 'dt':dt, \
        'cint':c0, 'D':D, 'alpha':alpha, 'kappa':kappa, 'K':K, 'delta':delta}

    #-- calculate total initial solute
    cint = np.sum(c0) 

    #-- solve numerically using matrix solver
    nsteps = 2

    for t in range(nsteps):
        c = rad_solver_2Daxi(c0, w, n, uw, vw, param_dict)
        #-- reformat solution to 2D
        c2d = np.zeros((ny, ny))
        for i in range(ny):
            for j in range(ny): 
                c2d[i, j] =  c[j*ny + i]  
        c0[:] = c2d[:]

    #-- calculate total final solute

    cfin = np.sum(c) 

    #-- check conservation

    assert cfin == approx(cint, rel = 1e-3)

    return