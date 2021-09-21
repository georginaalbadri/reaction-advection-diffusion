from Solver1D import PSOrad_solver_1D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from tqdm import tqdm
import pyswarms as ps
from Preliminary1D import InputVariables
from Preliminary1D import InputSoluteParameters
from PSOMinimiser1D import minimiser


#----------------------- SETUP MODEL PROBLEM 

#-- set problem geometry 
L = 5e-3  #length-scale [m]
T = 1e3  #timescale [s]

nx = 400 #number of grid points
dx = 1/nx #dimensionless grid spacing
dx2 = dx*dx

geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

#-- set variables
n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

#-- set initial condition and solute parameters
c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)


#-- combine variables and parameters in single dictionary for pyswarms function minimiser 
maindict = {'n':n, 'w':w, 'uw':uw, 'c0':c0, **param_dict}

#-- set time to run model
Tmax = 0.2
dt = param_dict['dt']
#add number of timesteps as parameter in dictionary 
nsteps = int(Tmax/dt)
maindict['nsteps'] = nsteps
print(nsteps)

# ---------------- NONDIMENSIONALISE PARAMETER BOUNDS (use as applicable)
#TODO: simplify this section 

# recover scalings/values needed
T = maindict['T']
L = maindict['L']
cM = 1e-9 
D = maindict['D']
D1 = D * L * L / T #recover dimensional D

# change dimensional minval, max val for each parameter as required to optimise within these bounds

#for D
minval, maxval = 1e-14, 1e-10 #minimum and maximum dimensional value
Dmin, Dmax = (T * minval) / (L**2), (T * maxval) / (L**2) #nondimensionalised

#for alpha
minval, maxval = 1e-13, 1e-11
almin, almax = (minval * L * L) / (D1 * cM), (maxval * L * L) / (D1 * cM)

#for kappa
minval, maxval = 1e-14, 1e-12 
kamin, kamax = (minval * L * L) / (D1 * cM), (maxval * L * L) / (D1 * cM)

#for delta
minval, maxval = 1e-5, 3e-4
delmin, delmax = minval * L * L / D1, maxval * L * L / D1

#for K
Kmin, Kmax = minval / cM, maxval / cM



#-------------------- CONFIGURE OPTIMISER

#-- select parameter(s) to optimise, and specify bounds 

#x = alpha, kappa - check this matches in PSOrad_solver_1D

dim = 2 #number of parameters to optimise
bounds = ([almin, almax], [kamin, kamax]) # (dimensionless) bounds of parameters to optimise
options = {'c1': 1.5, 'c2': 1.5, 'w': 0.5} #pyswarms global parameters 

# setup optimiser with above options 
optimizer = ps.single.GlobalBestPSO(n_particles = 20, dimensions = dim, options=options, bounds=bounds)


#----------------- RUN OPTIMISER - runs 'minimiser' function  

cost, pos = optimizer.optimize(minimiser, iters=100, n_processes = 4, **maindict)

print('cost = ', cost)
print('pos = ', pos)

with open('pywsarms_output.txt', 'a') as f:
    print(f'alpha, kappa = {pos}', file=f)
    print(f'cost = {cost}', file=f)




# ------------- Test optimal parameter values in function ------------------#

x = pos

#-- set problem geometry 
L = 5e-3  #length-scale [m]
T = 1e3  #timescale [s]

nx = 400 #number of grid points
dx = 1/nx #dimensionless grid spacing
dx2 = dx*dx

geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

#-- set variables
n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

#-- set initial condition and solute parameters
c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

maindict = maindict = {'n':n, 'w':w, 'uw':uw, 'c0':c0, **param_dict}

T = 0.2
dt = param_dict['dt']
nsteps = int(T/dt)
maindict['nsteps'] = nsteps


for t in tqdm(range(nsteps)):
    maindict = PSOrad_solver_1D(x, **maindict)

c = maindict['c0']
print('total solute =', np.sum(c))

