# reaction-advection-diffusion
1D, 2D cartesian and 2D axisymmetric solvers for reaction-advection-diffusion PDE. Also includes applications: parameter sweep, parameter sensitivity analysis (SALib), parameter optimisation (PSO - pyswarms). 

## Dependencies

Python 3.7 

numpy==1.20.3
matplotlib==3.0.3
tqdm==4.31.1
pytest==4.3.1
scipy==1.2.1
SALib==1.3.13

## Introduction

This repository contains solvers for a reaction-advection-diffusion PDE in 1D, 2D cartesian, and 2D axisymmetric (r-z). The solver files

- Solver1D.py

- Solver2D.py

- Solver2Daxi.py

each contain required functions to solve the equations at each timestep as a matrix problem, where the equation has been discretised using finite differences and is solver numerically using spsolve.  

The manager files

- Manager1D.py
- Manager2D.py

use the appropriate solvers to solve the PDE for given parameters, time and geometry.  

The parameter sweep function files 

ParameterSweep1D.py
ParameterSweep2D.py 

enable the functions to be solved for a range of parameters and/or timepoints to visually compare. 

Additionally, the repository contains two applications: Sensitivity Analysis and Parameter Optimisation. 

The files for sensitivity analysis

SensitivityAnalysis1D.py
SensitivityAnlaysis2D.py

use the Python library SALib to generate parameter samples, uses multiprocessing Pool to run these samples through the appriate solver, and uses a SALib analyser to generate a graph showing the sensitivity of a particular model outcome to each input parameter.

The files for parameter optimisation

PSO1D.py
PSO2D.py

use the Python library Pyswarms to optimise input parameters to match a set model outcome. 

[] $\alpha$

