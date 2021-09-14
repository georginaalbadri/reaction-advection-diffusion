import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import sys

from Solver1D import rad_solver_1D


#------------------------------------------#

# PARAMETER SWEEP

#------------------------------------------#
"""
    Parameters (input dimensional):
    c_int - initial uniform solute concentration [ng/ml]
    D - solute diffusion coefficient [m^2/s]
    alpha - solute production rate by cells [1/s]
    kappa - solution uptake rate by cells (Michaelis-Menten kinetics) [1/s]
    K - concentration of solute at which uptake is half-maximal [ng/ml]
    delta - solute decay rate (independent of cells) [1/s] 
"""

def ParamSweep1D(c0, n, w, uw, param_list, param_dict, param_name = 'D', Tmax = 1):

    dt = param_dict['dt']
    nsteps = int(Tmax / dt)
    print(nsteps)

    clist = []

    #-- run 1D model for each parameter value

    for val in tqdm(param_list):

        #change parameter value in dictionary to current sweep value (nondimensionalise)

        L = param_dict['L']
        T = param_dict['T']
        cM = 1e-9 #typical solute concentration [ng/ml]
    
        if param_name in ['D']:
            D = param_dict['D'] #extract original diffusion array
            D = (T * val) / (L**2) #replace diffusion in gel with new (non-dimensionalised) sweep value
            param_dict['D'] = D #reinstate new diffusion array

        if param_name in ['alpha', 'kappa']:
            D = param_dict['D']
            D1 = D * L * L / T
            param_dict[param_name] = (val * L * L) / (D1 * cM)

        if param_name in ['delta']:
            D = param_dict['D']
            D1 = D * L * L / T
            param_dict['delta'] = val * L * L / D1

        if param_name in ['K']:
            param_dict['K'] = val / cM

        #reset initial solute concentration
        c0 = param_dict['c0']

        for t in range(nsteps+1):
            c = rad_solver_1D(c0, n, w, uw, param_dict)

        #save result for each parameter value
        clist.append(c[:])


    #-- plot graph of results 
    
    nx = param_dict['nx']
    Y = np.linspace(0, 1, nx)

    plt.figure()

    # make legend handles
    num = len(param_list)
    legend = []
    for i in range(num):
        if np.abs(param_list[i]) > 99 or np.abs(param_list[i]) < 0.01:
            legend.append(f'${param_name}$ = {param_list[i]:.1E}')
        else:
            legend.append(f'${param_name}$ = {param_list[i]}')

    # make plot 
    plt.figure()
    for j in range(num):
        plt.plot(Y, clist[0+j], alpha = 0.7, lw = 1.5)
 
    if np.min(clist) >= 0:
        plt.ylim((0.9 * np.min(clist)), (1.1 * np.max(clist)))
    if np.min(clist) < 0:
        plt.ylim((1.1 * np.min(clist)), (1.1 * np.max(clist)))
        
    plt.xlabel('Height [-]')
    plt.ylabel('Solute concentration [-]')
    plt.legend(legend)
    plt.savefig(f'soluteconc_{param_name}_spatial.png', dpi = 400)

    return



