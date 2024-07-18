# Vitis HLS Source Code

This folder contains all the source files of the Vitis HLS project. To create the project, first download the files from this repository and copy them to a new folder. Open a terminal in that same folder and run the following command: “vitis_hls -f script.tcl”. After this, you can open Vitis HLS and simulate, using the provided testbench file, or synthesize the module to HDL code.

# Files contained:
- **l1_control.h:** This header file wraps all the required functions, variables, constants and libraries needed for the HoC module to work. All the constants are stored here, along with multiple pre-calculated divisions.

- **l1_control_sol1.cpp:** This file contains the HoC module. After receiving the current state as input data, it computes the sparsest control law and returns it as output data.

- **chol_decomp.h:** This file contains the Cholesky decomposition routine. This is used to compute the Cholesky factor corresponding to the KKT Matrix. 

- **ltris.h:** This file contains the LTRIS (Lower Triangular Solver) and UTRIS (Upper Triangular Solver) routines, used for solving the triangular systems of equations resulted after the Cholesky factorization.

- **hoc_sim_tb.cpp:** This file is the testbench of the project. It can also be used as a simulator. After providing the current state of the system to the control algorithm, it receives the control signals, applies them, simulates the trajectory of the system and repeats the process for a predefined number of steps. The simulation data, stored in files “state_file.dat” and “hardware_output_file.dat”, can be further copied and then represented using the MATLAB script “parse_plot.m”.


