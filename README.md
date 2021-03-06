[![Build Status](https://app.travis-ci.com/georginaalbadri/reaction-advection-diffusion.svg?branch=main)](https://app.travis-ci.com/georginaalbadri/reaction-advection-diffusion)

# reaction-advection-diffusion

Author: Georgina Al-Badri
georgina.kennedy.13@ucl.ac.uk

1D and 2D axisymmetric solvers for reaction-advection-diffusion PDE. Also includes applications: parameter sweep, parameter sensitivity analysis (SALib), parameter optimisation (PSO - pyswarms). 

## Introduction

### ExampleCase.ipynb (1D)

This Jupyter notebook displays outputs for an example use of the following code/applications contained in this repository. 

### Solvers

This repository contains solvers for a reaction-advection-diffusion PDE in 1D and 2D axisymmetric (r-z axis). The solver files

- Solver1D.py
- Solver2D.py (axisymmetric)

each contain required functions to solve the equations at each timestep as a matrix problem, where the equation has been discretised using finite differences and is solver numerically using spsolve.  

Additionally Solver1D.py contains the solver PSOrad_solver_1D, which is equivalent to rad_solver_1D except that it takes only a dictionary argument, for use with particle swarm optimisation (called in PSOMinimiser1D and PSO1D). 

The parameter sweep function file

- ParameterSweep1D.py

is an alternative for the 1D model, based on Solver1D, to enable the function to be solved for a range of parameters and/or timepoints to visually compare outcomes. 

### Managers

The manager files

- Manager1D.py
- Manager2D.py

use the appropriate solvers to solve the PDE for given parameters, time and geometry.  


### Applications 

Additionally, the repository contains two applications: Sensitivity Analysis and Parameter Optimisation. 

The file for sensitivity analysis

- SensitivityAnalysis1D.py

uses the Python library SALib to generate parameter samples, uses multiprocessing Pool to run these samples through the appriate solver, and uses a SALib analyser to generate a graph showing the sensitivity of a particular model outcome to each input parameter.

The files for parameter optimisation

- PSO1D.py
- PSOMinimiser1D.py

use the Python library Pyswarms to optimise input parameters to match a set model outcome. The PSO is run through PSO1D, which minimised the function in Minimiser1D (based on Solver1D.py function PSOrad_solver_1D.py). 

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

Manager*.py 

Requires Solver and Preliminary 

### Supporting files

Solver*.py

Parameter Sweep*.py (uses Solver*.py)

### Application execution files

SensitivityAnalysis*.py and PSO*.py

### Test files

## Updates

September 2021: All solvers and applications in 1D ready and tested. Next update will bring 2D file structure in line with 1D, and add sesntivity analysis and parameter sweep applications.

Addition of flux boundary conditions projected October 2021. 

