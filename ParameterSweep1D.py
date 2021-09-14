import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import sys

from Preliminary1D import InputVariables
from Preliminary1D import InputSoluteParameters
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

def ParamSweep1D(c0, n, w, uw, param_list, param_dict, tlist, param_name = 'chi', Tmax = 1, plot_type = 'space'):

    dt = param_dict['dt']
    print(dt)
    nsteps = int(Tmax / dt)

    print(nsteps)

    clist = []
    ctot = []

    steplist = []
    for i in range(len(tlist)):
        steplist.append(int(tlist[i] / dt))

    #-- Run 1D model for each parameter value

    for val in tqdm(param_list):

        #change parameter value in dictionary to current sweep value (first nondimensionalise)

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

        ctime = []

        #reinitialize variables
        c0 = param_dict['c0']

        for t in range(nsteps+1):
            c = rad_solver_1D(c0, n, w, uw, param_dict)

            nx = param_dict['nx']
            
            ctime.append(np.sum(c[:nx])/nx)

            if t in steplist:
                clist.append(c[:])

        ctot.append(ctime)


    #-- Plot graph of variables

    if plot_type in ['space']:
    
        nx = param_dict['nx']
        Y = np.linspace(0, 1, nx)

        plt.figure()

        num = len(param_list)

        legend = []
        for i in range(num):
            if np.abs(param_list[i]) > 99 or np.abs(param_list[i]) < 0.01:
                legend.append(f'$\{param_name}$ = {param_list[i]:.1E}')
            else:
                legend.append(f'$\{param_name}$ = {param_list[i]}')

        q = len(tlist)
        if q>4 or q<4:
            sys.exit('number of timepoints must be 4 to do suplot')

        plt.figure()
        fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)

        for j in range(num):
            axs[0, 0].plot(Y, clist[0+j*q], alpha = 0.7, lw = 1.5)
            axs[0, 0].set_title(f'T = {tlist[0]}')
            axs[0, 1].plot(Y, clist[1+j*q], alpha = 0.7, lw = 1.5)
            axs[0, 1].set_title(f'T = {tlist[1]}')
            axs[1, 0].plot(Y, clist[2+j*q], alpha = 0.7, lw = 1.5)
            axs[1, 0].set_title(f'T = {tlist[2]}')
            axs[1, 1].plot(Y, clist[3+j*q], alpha = 0.7, lw = 1.5)
            axs[1, 1].set_title(f'T = {tlist[3]}')
            if np.min(clist) >= 0:
                axs[0, 0].set_ylim((0.9 * np.min(clist)), (1.1 * np.max(clist)))
            if np.min(clist) < 0:
                axs[0, 0].set_ylim((1.1 * np.min(clist)), (1.1 * np.max(clist)))
            for ax in axs.flat:
                ax.set(xlabel='Height [-]', ylabel= 'Solute concentration [-]')
            # Hide x labels and tick labels for top plots and y ticks for right plots.
            for ax in axs.flat:
                ax.label_outer()

        plt.subplots_adjust(right=0.75)
        fig.legend(legend, bbox_to_anchor = (0.98, 0.9))

        #bbox_to_anchor = (1, 0.9)
        
        plt.savefig(f'soluteconc_{param_name}_subplot.png', dpi = 400)
        #bbox_inches='tight'

    if plot_type in ['time']:

        dt = param_dict['dt']
        nsteps = int(Tmax / dt)
        T1 = np.linspace(dt, Tmax, nsteps)

        # Number of parameter values tested
        num = len(ctot)

        # Legend 
        legend = []
        for i in range(num):
            if np.abs(param_list[i]) > 99 or np.abs(param_list[i]) < 0.01:
                legend.append(f'${param_name}$ = {param_list[i]:.1E}')
            else:
                legend.append(f'${param_name}$ = {param_list[i]}')

    # Plot time profile for each param value on same figure

        plt.figure()

        for i in range(num):
            plt.plot(T1, ctot[i][1:], alpha = 0.7, lw = 1)
            plt.xlabel('Time [-]')
            plt.ylabel('Solute concentration [-]')
            plt.title('Average solute concentration over time')
            plt.legend(legend)

        plt.savefig(f'soluteconc_{param_name}_timeplot.png', dpi = 400, bbox_inches = 'tight')

    return



