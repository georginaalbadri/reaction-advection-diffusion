import numpy as np
import random
import sys



#---------------------------------------------------#

# CELL/WATER VOLUME DISTRIBUTIONS

#---------------------------------------------------#


def InputVariables(parameters_dict, n_option = "random", nmin = 0.1, nmax = 0.2, m = 0.03):
    """
    Returns distributions for cells, water (solvent), and water (solvent) velocity 

    About 'n_option': Input takes 'uniform', 'random', 'linear', or 'sinusoidal'. 
    - Uniform gives initial cell volume fraction 'nmin' across whole domain.
    - Random produces array of values 'nmin' to 'nmax' randomly distributed.
    - Linear produces gradient distribution from 'nmin' at LHS to 'nmax' at RHS
    - Sinusoidal produces an oscillating distribution with minima at 'nmin' and maxima at 'nmax' 

    Variables:
    n - cell volume fraction [-]
    w - water volume fraction [-] 
    uw - water velocity [m/s] 
    NB: next code update will include advection of water in and out of the domain 

    Parameters:
    nmin, nmax - minimum and maximum cell volume fraction in distribution, used to generate cell distribution depending on 'n_option'
    m - solid matrix volume fraction, between 0 and 1 [-]
    phi - available volume fraction in gel, phi = 1-m [-]

    """
    nx = parameters_dict['nx']
    dx = parameters_dict['dx']

    # set cell initial distribution based on function input
    while n_option not in ['uniform', 'random', 'linear', 'sinusoidal']:
        print("Invalid initial cell distribution choice made (can be 'uniform', 'random', 'linear' or 'sinusoidal')")
        exit()

    if n_option in ['uniform']: #selects uniform distribution n = nmin 
        n = nmin * np.ones((nx))

    if n_option in ['random']: #selects distribution with random fluctuations between cmin and cmax
        np.random.seed(42)
        n = nmin + ((nmax - nmin) * np.random.rand(nx))
    
    if n_option in ['linear']: #selects linear distribution between cmin and cmax
        n = np.zeros((nx))
        for i in range(nx):
            n[i] = nmin + ((nmax - nmin) / (nx-1)) * (i)
    
    if n_option in ['sinusoidal']:
        n = (nmin + ((nmax - nmin) / 2)) * np.ones((nx))
        for i in range(nx):
            n[i] += ((nmax - nmin) / 2) * np.sin(20 * np.pi * i * dx)

    # amount of free volume
    phi = 1 - m

    # water volume fraction dependent on cell distribution via no voids constraint (n + w + m = 1)
    w = phi - n 

    # water velocity
    uw = np.zeros((nx))

    parameters_dict['phi'] = phi
    parameters_dict['m'] = m

    return n, w, uw, parameters_dict


#------------------------------------------------------#

# SOLUTE INITIAL CONDITIONS AND PARAMETERS

#------------------------------------------------------#
    

def InputSoluteParameters(parameters_dict, c_int = 0, D = 1e-10, alpha = 1e-12, kappa = 1e-13, K = 0.1, delta = 2.5*1e-4, dt_mult = 10):
    
    """
    Updates variable dictionary with solute concentration initial conditions, and updates parameter dictionary with solute parameter values

    Variables and initial conditions:
    c0 - initial solute concentration [ng/ml]

    Parameters (input dimensional):
    c_int - initial uniform solute concentration [ng/ml]
    D - solute diffusion coefficient [m^2/s]
    alpha - solute production rate by cells [1/s]
    kappa - solution uptake rate by cells (Michaelis-Menten kinetics) [1/s]
    K - concentration of solute at which uptake is half-maximal [ng/ml]
    delta - solute decay rate (independent of cells) [1/s] 

    """
    nx = parameters_dict["nx"] #retrieve grid size 
    dx = parameters_dict['dx']

    #-- set initial solute concentration based on cint choice 
    c0 = c_int * np.ones((nx))

    #-- nondimensionalise parameters
    L = parameters_dict['L']
    T = parameters_dict['T']
    cM = 1e-9 #concentration ng/ml

    D1 = D * T / L**2

    alpha1 = (alpha * T / cM) #non-dim production rate 
    kappa1 = (kappa * T / cM) #non-dim production rate 
    delta1 = (delta * T) #non-dim degradation rate 

    #-- set appropriate timestep based on explicit FD scheme limits 
    #dt_mult scales this timestep 
    dx2 = parameters_dict['dx2']
    maxtimestep = (dx2 / (2 * D1) )
    dt = dt_mult * maxtimestep

    #-- update parameters dictionary
    parameters_dict["dt"] = dt
    parameters_dict['c0'] = c0 
    parameters_dict["D"] = D1
    parameters_dict["alpha"] = alpha1
    parameters_dict["kappa"] = kappa1
    parameters_dict["K"] = K
    parameters_dict["delta"] = delta1

    return c0, parameters_dict 
