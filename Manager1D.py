
import numpy as np
import matplotlib.pyplot as plt

from Preliminary1D import InputVariables
from Preliminary1D import InputSoluteParameters
from Solver1D import rad_solver_1D

import time 
from tqdm import tqdm



#------------------------------------------#

# Run 1D reaction-advection-diffusion model

#------------------------------------------#

#-- set problem geometry 
L = 5e-3  #length-scale [m]
T = 1e5  #timescale [s]

nx = 500 #number of grid points
dx = 1/nx #dimensionless grid spacing
dx2 = dx*dx

geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

#-- set variables
n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "random", nmin = 0.1, nmax = 0.2, m = 0.03)

#-- set initial condition and solute parameters
c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 1, D = 1e-10, dt_mult = 100)


#-- plot initial distribution
#plt.figure()
#plt.plot(c0)
#plt.title('Initial solute distribution')
#plt.savefig('initial_solute_distribution.png', dpi = 300)


#-- calculate number of timesteps to run 
Tmax = 0.1 #maximum time to run [1e5 seconds]
dt = param_dict['dt']
nsteps = int(Tmax/dt)

print('no. of timesteps', nsteps)

#-- track time taken to run 
start_time = time.time() 


#-- iterate RAD solver over number of timesteps

for t in tqdm(range(nsteps)):
    c = rad_solver_1D(c0, n, w, uw, param_dict)
    c0 = c.copy() #update previous timestep value


#-- plot initial and final distribution on same graph
cint = param_dict['c0']

plt.figure()
plt.plot(cint, 'b')
plt.plot(c, 'g')
plt.legend(('T = 0', f'T = {Tmax}'))
plt.title(f'Solute distribution at T = 0 and T = {Tmax}')
plt.savefig(f'solute_dist_T{Tmax}.png', dpi = 300)


"""
#--------------------------------------------------------------------#

# Run parameter sweep for 1D reaction-advection-diffusion model

#--------------------------------------------------------------------#

#-- set problem geometry 
L = 5  #length-scale [mm]
T = 1e5  #timescale [s]

nx = 500 #number of grid points
dx = 1/nx #dimensionless grid spacing
dx2 = dx*dx

geometry_dict =  {'L':L, 'T':T, 'nx':nx, 'dx':dx, 'dx2':dx2} 

#-- set variables
n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "random", nmin = 0.1, nmax = 0.2, m = 0.03)

#-- set initial condition and solute parameters
c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-11, alpha = 1e-12, kappa = 1e-13, K = 0.1, delta = 2.5*1e-4, dt_mult = 1000)



#-- run parameter sweep for specified parameter
ParamSweep(paramlist, var_dict, param_dict, tlist, param_name = 'D1', T = 10, plot_type = 'time, space', subplot = not None, save = not None)
"""