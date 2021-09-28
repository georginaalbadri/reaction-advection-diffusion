
import numpy as np
import matplotlib.pyplot as plt
import time 
from tqdm import tqdm
from Preliminary1D import InputGeometry, InputVariables
from Preliminary1D import InputSoluteParameters
from Solver1D import rad_solver_1D, rad_solver_flow_1D
from ParameterSweep1D import ParamSweep1D, ParamSweepPlot1D


#-------------------------------------------------#

# Run 1D reaction-advection-diffusion model

#-------------------------------------------------#

geometry_dict = InputGeometry(L = 1e-3, T = 1e3, nx = 400) 

n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)


#-- plot initial distribution
plt.figure()
plt.plot(c0)
plt.title('Initial solute distribution')
plt.savefig('initial_solute_distribution.png', dpi = 300)


#-- calculate number of timesteps to run 
Tmax = 0.5 #maximum (dimensionless) time to run 
dt = param_dict['dt']
nsteps = int(Tmax/dt)

print('no. of timesteps', nsteps)

#-- track time taken to run 
start_time = time.time() 


#-- iterate RAD solver over number of timesteps
for t in tqdm(range(nsteps)):
    c = rad_solver_flow_1D(c0, n, w, uw, param_dict)
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

geometry_dict = InputGeometry(L = 1e-3, T = 1e3, nx = 400) 

n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = 1e-10, dt_mult = 10)


param_list = [1e-13, 1e-12, 1e-11, 5e-11]

Tmax = 1
dt = param_dict['dt']
nsteps = int(Tmax / dt)

#-- run parameter sweep for specified parameter values
clist = ParamSweep1D(c0, n, w, uw, param_list, param_dict, param_name = 'alpha', nsteps = nsteps)

#-- plot results on single figure 
ParamSweepPlot1D(param_dict, clist, param_list, param_name = 'alpha')
"""