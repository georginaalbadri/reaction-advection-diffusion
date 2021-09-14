import numpy as np
import random
import sys



#---------------------------------------------------#

# CELL/WATER VOLUME DISTRIBUTIONS

#---------------------------------------------------#


def InputVariables(parameters_dict, n_option = "random", nmin = 0.1, nmax = 0.2, m = 0.03):
    """
    Create arrays for other species (cell distribution, water distribution, and water velocity) 
    Returns dictionary of variables, and updates dictionary of parameters 

    About 'n_option': Input takes 'uniform', 'random', 'linear', or 'sinusoidal'. 
    - Uniform gives initial cell volume fraction 'nmin' across whole domain.
    - Random produces array of values 'nmin' to 'nmax' randomly distributed.
    - Linear produces gradient distribution from 'nmin' at LHS to 'nmax' at RHS
    - Sinusoidal produces an oscillating distribution with minima at 'nmin' and maxima at 'nmax' 

    Variables:
    n - cell volume fraction [-]
    w - water volume fraction [-] 
    uw, vw - x, y components of water velocity [m/s]

    Parameters:
    nmin, nmax - minimum and maximum cell volume fraction in distribution, used to generate cell distribution depending on 'n_option'
    m - solid matrix volume fraction, between 0 and 1 [-]
    phi - available volume fraction in gel, phi = 1-m [-]

    """
    nx, ny = parameters_dict['nx'], parameters_dict['ny'] #retrieve grid size
    dx = parameters_dict['dx']

    # set cell initial distribution based on function input
    while n_option not in ['uniform', 'random', 'linear', 'sinusoidal']:
        print("Invalid initial cell distribution choice made (can be 'uniform', 'random', 'linear' or 'sinusoidal')")
        exit()

    if n_option in ['uniform']: #selects uniform distribution n = nmin 
        n = nmin * np.ones((nx, ny))

    if n_option in ['random']: #selects distribution with random fluctuations between cmin and cmax
        np.random.seed(42)
        n = nmin + ((nmax - nmin) * np.random.rand(nx, ny))
    
    if n_option in ['linear']: #selects linear distribution between cmin and cmax
        n = np.zeros((nx, ny))
        for i in range(ny):
            n[i, :] = nmin + ((nmax - nmin) / (ny-1)) * (i)
    
    if n_option in ['sinusoidal']:
        n = (nmin + ((nmax - nmin) / 2)) * np.ones((nx, ny))
        for i in range(ny):
            n[i, :] += ((nmax - nmin) / 2) * np.sin(20 * np.pi * i * dx)

    # amount of free volume
    phi = 1 - m

    # water volume fraction dependent on cell distribution via no voids constraint (n + w + m = 1)
    w = phi - n 

    # water velocity 
    uw = np.zeros((nx, ny))
    vw = np.zeros((nx, ny))

    # create variables dictionary

    # update parameters dictionary 
    parameters_dict["phi"] = phi
    parameters_dict["m"] = m 

    return n, w, uw, vw, parameters_dict


#----------------------------------#

# VEGF PRELIMINARIES

#----------------------------------#
    

def InputSoluteParameters(parameters_dict, c_int = 0, D = 1e-11, alpha = 1e-11, kappa = 1e-13, K = 0.1, delta = 2.5*1e-4, dt_mult = 1000, dt = 0.001):
    
    """
    Updates variable dictionary with solute concentration initial conditions, and updates parameter dictionary with solute parameter values

    Variables and initial conditions:
    c - solute concentration [ng/ml]
    c0 - initial solute concentration in gel [ng/ml]

    Parameters (input dimensional):
    cint = initial uniform concentration [ng/ml]
    D - solute diffusion coefficient [m^2/s]
    alpha - solute production rate by cells [1/s]
    kappa - solution uptake rate by cells (Michaelis-Menten kinetics) [1/s]
    K - concentration of solute at which uptake is half-maximal [ng/ml]
    delta - solute decay rate (independent of cells) [1/s] 

    """
    nx, ny = parameters_dict["nx"], parameters_dict["ny"] #retrieve grid size 
    dx = parameters_dict['dx']

    #-- set initial solute concentration based on c_option choice 
    c0 = c_int * np.ones((nx, ny))

    #-- nondimensionalise parameters
    L = parameters_dict['L']
    T = parameters_dict['T']
    cM = 1e-9 #concentration ng/ml

    Pe = L**2 / (T * D) #non-dim Peclet number

    alpha1 = (alpha * L**2) / (D * cM) #non-dim production rate 
    kappa1 = (kappa * L**2) / (D * cM) #non-dim production rate 
    delta1 = (delta * L**2) / D #non-dim degradation rate 

    #-- set appropriate timestep based on explicit FD scheme limits 
    #dt_mult scales this timestep 
    dx2 = parameters_dict['dx2']
    maxtimestep = (dx2 * Pe / 2 )
    dt = dt_mult * maxtimestep
    dt = dt

    #-- update parameters dictionary
    parameters_dict["dt"] = dt
    parameters_dict['c0'] = c0 
    parameters_dict["Pe"] = Pe
    parameters_dict["alpha"] = alpha1
    parameters_dict["kappa"] = kappa1
    parameters_dict["K"] = K
    parameters_dict["delta"] = delta1

    return c0, parameters_dict 


