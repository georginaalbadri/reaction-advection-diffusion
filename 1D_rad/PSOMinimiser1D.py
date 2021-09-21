import numpy as np
from scipy.ndimage import label
from tqdm import tqdm
import pyswarms as ps
from Solver1D import PSOrad_solver_1D


@ps.cost
def minimiser(x, **maindict):

    #-- run model 
    nsteps = maindict['nsteps']
    for t in range(nsteps):
        maindict = PSOrad_solver_1D(x, **maindict)

    #-- compare model result to target value, to minimise difference 

    c = maindict['c0']

    # total solute concentration
    ctot = np.sum(c)

    # target solution concentration
    creal = 7.5

    # minimise difference between real and simulated total concentration
    minimise = abs(ctot - creal)

    return minimise 
    
    