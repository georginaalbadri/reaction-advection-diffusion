from Solver1D import PSOrad_solver_1D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from tqdm import tqdm
import pyswarms as ps
from Preliminary1D import InputGeometry, InputVariables, Nondimensionalise, InputSoluteParameters
from PSOFunctions1D import minimiser


#----------------------- SETUP MODEL PROBLEM 

geometry_dict = InputGeometry(L = 1e-3, T = 1e3, nx = 400) 

n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)


# combine variables and parameters in single dictionary for pyswarms function
maindict = {'n':n, 'w':w, 'uw':uw, 'c0':c0, **param_dict}

#-- set time to run model
Tmax = 0.2
dt = param_dict['dt']
#add number of timesteps as parameter in dictionary 
nsteps = int(Tmax/dt)
maindict['nsteps'] = nsteps
print(nsteps)

# ---------------- NONDIMENSIONALISE PARAMETER BOUNDS (change parameter name/bounds as applicable)

parameters = ['alpha', 'kappa']
dim_bounds = [[1e-13, 1e-11], [1e-14, 1e-12]]

min1 = Nondimensionalise(parameter_dict, dimval = dim_bounds[0][0], param_name = parameters[0])
max1 = Nondimensionalise(parameter_dict, dimval = dim_bounds[0][1], param_name = parameters[0])
min2 = Nondimensionalise(parameter_dict, dimval = dim_bounds[1][0], param_name = parameters[1])
max2 = Nondimensionalise(parameter_dict, dimval = dim_bounds[1][1], param_name = parameters[1])


#-------------------- CONFIGURE OPTIMISER

dim = len(parameters)
bounds = ([min1, max1], [min2, max2]) #(dimensionless) bounds from above
options = {'c1': 1.5, 'c2': 1.5, 'w': 0.5} #pyswarms global parameters 

optimizer = ps.single.GlobalBestPSO(n_particles = 20, dimensions = dim, options=options, bounds=bounds)


#----------------- RUN OPTIMISER - runs 'minimiser' function in particle swarm 

cost, pos = optimizer.optimize(minimiser, iters=100, n_processes = 4, **maindict)

print('cost = ', cost)
print('pos = ', pos)

with open('pywsarms_output.txt', 'a') as f:
    print(f'alpha, kappa = {pos}', file=f)
    print(f'cost = {cost}', file=f)




# -------------------- TEST OUTPUT -------------------------#

x = pos

geometry_dict = InputGeometry(L = 1e-3, T = 1e3, nx = 400) 

n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)

maindict = maindict = {'n':n, 'w':w, 'uw':uw, 'c0':c0, **param_dict}

nsteps = int(Tmax/dt)
maindict['nsteps'] = nsteps

for t in tqdm(range(nsteps)):
    maindict = PSOrad_solver_1D(x, **maindict)

c = maindict['c0']
print('total solute =', np.sum(c))

