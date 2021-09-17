# reaction-advection-diffusion

Georgina Al-Badri
zcahgke@ucl.ac.uk


1D, 2D cartesian and 2D axisymmetric solvers for reaction-advection-diffusion PDE. Also includes applications: parameter sweep, parameter sensitivity analysis (SALib), parameter optimisation (PSO - pyswarms). 

# Model

$\frac{\partial c w }{\partial t}$

## Dependencies

Python 3.7 

numpy==1.20.3
matplotlib==3.0.3
tqdm==4.31.1
pytest==4.3.1
scipy==1.2.1
SALib==1.3.13

## Introduction

### Solvers

This repository contains solvers for a reaction-advection-diffusion PDE in 1D, 2D cartesian, and 2D axisymmetric (r-z). The solver files

- Solver1D.py
- Solver2D.py
- Solver2Daxi.py

each contain required functions to solve the equations at each timestep as a matrix problem, where the equation has been discretised using finite differences and is solver numerically using spsolve.  

The parameter sweep function files 

- ParameterSweep1D.py
- ParameterSweep2D.py 

are an alternative, based on the Solvers, to enable the functions to be solved for a range of parameters and/or timepoints to visually compare outcomes. 

### Managers

The manager files

- Manager1D.py
- Manager2D.py

use the appropriate solvers to solve the PDE for given parameters, time and geometry.  


### Applications 

Additionally, the repository contains two applications: Sensitivity Analysis and Parameter Optimisation. 

The files for sensitivity analysis

- SensitivityAnalysis1D.py
- SensitivityAnlaysis2D.py

use the Python library SALib to generate parameter samples, uses multiprocessing Pool to run these samples through the appriate solver, and uses a SALib analyser to generate a graph showing the sensitivity of a particular model outcome to each input parameter.

The files for parameter optimisation

- PSO1D.py
- PSO2D.py

use the Python library Pyswarms to optimise input parameters to match a set model outcome. 

### Variables

- `c` : solute 
- `u`, `v` : x, y components of solvent velocity (1D/2D cartesian)
- `u`, `v` : z, r components of solvent velocity (2D axisymmetric)
- `w` : water volume fraction (w=1 if modelling solute in liquid; w<1 if modelling solute in liquid within solid/rigid porous scaffold). 

### Parameters
- `D` : diffusion coefficient of solute
- `alpha` : production rate of solute
- `delta` : degradation rate of solute
- `kappa` : uptake rate of solution (Michaelis-Menten kinetics)
- `K` : concentration of solute at which uptake is half-maximal (Michaelis-Menten kinetics)
- `dt` : timestep
- `h` : grid spacing 
- `L` : height/length of model geometry 
- `T` : typical time 

## File structure

### Execution files

Manager*.py OR Parameter Sweep*.py

Requires Solver and Preliminary 

### Supporting files

Solver*.py

### Application execution files

SensitivityAnalysis*.py and PSO*.py



