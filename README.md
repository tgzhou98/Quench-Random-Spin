# Quench-Random-Spin

## Contents

This repository contains codes for the paper *Emergent Universal Quench Dynamics in Randomly Interacting Spin Models*. There are three packages for theoretical simulations and one package for producing the figure.

- Packages

  - **largeM_dynamics**: Numerical simulation of the large-$M$ dynamics.
  
  - **exact_diagonalization**: Numerical simulation of the exact-diagonalization dynamics.
  
  - **meanfield**: Numerical simulation of the mean-field dynamics.
  
  - **MATLAB_fit**: Serious analysis of error bar and best fitting parameter.

## Dependence

- Julia
  - Julia 1.7.1 (or above)
  - Packages: `LinearAlgebra`, `JLD`, `JLD2`, `DelimitedFiles`, `Plots`, `Distributed`, `Distributed`, `DataFrames`, `MAT`, `Random`, `Statistics`, `ProgressMeter`
  - Simulate Large-M, exact-diagonalization and mean-field numerics.

- Mathematica 12.1 (or above)

- Matlab R2022b (or above)
  - Curve Fitting Toolbox

## Running the program

- **largeM_dynamics**
    - Run `job_script_dynamics.jl` to simulate the large-$M$ dynamics. Change the parameters in the file if you want to investigate other cases.
    - Run `analyze_dynamics_anyf_for_exp_draft.nb` to plot the magnetization evolution and simply fit the curve. It also exports data to the Matlab.
    - `dynamics_20230715_draftRandomf_hightemp_SpinModel_Jbar0_rhox` folder contains parts of the simulation result.

-  **mean_field**
    - Run `meanfield_randomspin.jl` to simulate the mean-field dynamics. Change the parameters in the file if you want to investigate other cases.
    - Run `meanfield_data_to_MATLAB.jl` to export the time evolution data to Matlab datatype.
  `figs_exp_x` folder contains the simulation result.

-  **exact_diagonalization**
    - Run `dynamics_quench_Sx.jl` to simulate the exact diagonalization dynamics. Change the parameters in the file if you want to investigate other cases. This program automatically exports to Matlab datatype.
    - `rs_Sx_draftRandomf_dynamics_20230716_rhox` folder contains the simulation result.

- **MATLAB_fit**
    - `Sx` folder contains data for three numerical simulations, and `Exp_PlotData_new_method` folder is for experimental data.
    - `LargeM_quench_x_all.mlx`, `Meanfield_quench_x_all.mlx`, `ED_quench_x_all.mlx`, `Exp_quench_x_all.mlx` are four separate Matlab notebook files to serious deal with the standard deviation and best fitting parameter.