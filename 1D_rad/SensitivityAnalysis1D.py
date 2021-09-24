import matplotlib.pyplot as plt
import numpy as np
import time
import sys
from tqdm import tqdm 
import multiprocessing as mp
from multiprocessing import Pool
from SALib.sample import saltelli
from SALib.analyze import sobol
from Preliminary1D import InputGeometry, InputVariables, InputSoluteParameters
from Solver1D import rad_solver_1D


"""
    Parameters :
    c_int - initial uniform solute concentration [ng/ml]
    D - solute diffusion coefficient [m^2/s]
    alpha - solute production rate by cells [1/s]
    kappa - solution uptake rate by cells (Michaelis-Menten kinetics) [1/s]
    K - concentration of solute at which uptake is half-maximal [ng/ml]
    delta - solute decay rate (independent of cells) [1/s] 

"""

# ---------------------------- SA SETUP -------------------------- #

# set parameter value bounds
problem = {
    'num_vars': 5,
    'names': ['D', 'alpha', 'kappa', 'K', 'delta'],
    'bounds': [[1e-12, 1e-9],
               [1e-13, 1e-11],
               [1e-15, 1e-12],
               [0.05, 0.2],
               [1e-5, 5e-4]]
}

# use saltelli sampler to generate parameter samples within bounds given 
# number argument in sample should be large as possible (by <1000) for high confidence values 
param_values = saltelli.sample(problem, 512)
nparams = len(param_values)
print(f'number of parameter samples to run = {nparams}')

# --------------------- MULTIPROCESSING FUNCTION ---------------------- #

def multiprocessing_func(param_values_sample):

    #unpack parameter values from sampler
    D = param_values_sample[0]
    alpha = param_values_sample[1]
    kappa = param_values_sample[2]
    K = param_values_sample[3]
    delta = param_values_sample[4]

    #-- setup model
    geometry_dict = InputGeometry(L = 5e-3, T=1e3, nx = 200)

    n, w, uw, parameter_dict = InputVariables(geometry_dict, n_option = "sinusoidal", nmin = 0.1, nmax = 0.2, m = 0.03)

    c0, param_dict = InputSoluteParameters(parameter_dict, c_int = 0, D = D, alpha = alpha, kappa = kappa, K = K, delta = delta, dt_mult = 10)

    #-- calculate number of timesteps
    dt = param_dict['dt']
    Tmax = 1
    nsteps = int(Tmax/dt)
    
    assert nsteps != 0 


    ##-- optional - to analyse multiple timepoints
    #tnum = 2 #number of timepoints to capture sensitivity 
    #tlist = np.linspace(1, nsteps-1, tnum)
    #tlist = [int(item) for item in tlist]
    #ctotlist = []

    #run timesteps
    for t in range(nsteps):
        c = rad_solver_1D(c0, n, w, uw, param_dict)
        c0 = c.copy()

        ## optional - conduct sensitivity analysis at multiple timepoints
        #if t in tlist:
            ##save total solute
            #ctotlist.append(sum(c)) 
    

    # for single output at Tmax
    ctot = np.sum(c[:])

    # for multiple timepoints, return ctotlist
    return ctot


# ------------------------ CODE EXECUTION -------------------------- #

if __name__ == '__main__':

    # choose number of processes to run based on number of cpus of machine 
    num_cpu = mp.cpu_count() - 1
    print(f'Running on {num_cpu} threads')

    start_time = time.time()

    #-- run model across num_cpu processes using multiprocessing.Pool module
    with Pool(processes = num_cpu) as pool:

        # input function to run, and list of parameter values param_values to run it for
        # returns list of model outputs (ctot) in order of samples passed 
        results = list(tqdm(pool.imap(multiprocessing_func, param_values), total = nparams))
   

    # determine time takem 
    print('time taken', time.time() - start_time)

    # convert list output to array for analysis by SALib module
    results_set = np.array(results)

    #-- for MULTIPLE TIMEPOINTS (ctotlist), analyse one time point at a time 

    ## transpose and save output
    #results_set = np.transpose(results_set)
    #np.save('results_set', results_set)

    # analyse data for each time point 
    #ST_set = []   
    #ST_conf = []   

    #for result in results_set:
    #    ST = sobol.analyze(problem, result)
    #    ST_set.append(ST['ST'])
    #    ST_conf.append(ST['ST_conf'])

    ## transpose back and save sensitivity values and confidence levels 
    #ST_set=np.array(ST_set)
    #ST_set=np.transpose(ST_set)
    #np.save('ST_set_array', ST_set)

    #ST_conf=np.array(ST_conf)
    #ST_conf=np.transpose(ST_conf)
    #np.save('ST_conf_array', ST_conf)

    #-- for SINGLE TIMEPOINT (ctot) output

    # analyse results_set, array of model outputs 

    ST = sobol.analyze(problem, results_set, print_to_console=False)

    ST_set = ST['ST']
    ST_conf = ST['ST_conf']

    #print list of parameter names
    print(problem['names'])
    #print list of sensitivity values (in order as above)
    print(ST_set)
    #print list of confidence values 
    print(ST_conf)

    #----------------- PLOT SA RESULTS - optional - for MULTIPLE TIMEPOINTS only  ----------------#
    #z = np.linspace(0, 1, 3)
#
    ## cell parameters plot 
    #fig, ax = plt.subplots()
#
    ## plot sensitivity values
    #ax.plot(z, ST_set[0], label = 'D')
    #ax.plot(z, ST_set[1], label = r'$\alpha$')
    #ax.plot(z, ST_set[2], label = r'$\kappa$')
    #ax.plot(z, ST_set[3], label = 'K')
    #ax.plot(z, ST_set[4], label = r'$\delta$')
#
    ## plot confidence intervals 
    #ax.fill_between(z, ST_set[0]- ST_conf[0], ST_set[0]+ ST_conf[0], alpha=0.2)
    #ax.fill_between(z, ST_set[1]- ST_conf[1], ST_set[1]+ ST_conf[1], alpha=0.2)
    #ax.fill_between(z, ST_set[2]- ST_conf[2], ST_set[2]+ ST_conf[2], alpha=0.2)
    #ax.fill_between(z, ST_set[3]- ST_conf[3], ST_set[3]+ ST_conf[3], alpha=0.2)
    #ax.fill_between(z, ST_set[4]- ST_conf[4], ST_set[4]+ ST_conf[4], alpha=0.2)
#
    #ax.legend()
    #plt.savefig('SA_rad_parameters_time.png')
    #plt.show()
#
#
#