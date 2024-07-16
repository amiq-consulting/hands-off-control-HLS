############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
## Copyright 2022-2023 Advanced Micro Devices, Inc. All Rights Reserved.
############################################################
open_project Hands_off_control_HLS
set_top hoc_module
add_files Hands_off_control_HLS/chol_decomp.h
add_files Hands_off_control_HLS/l1_control.h
add_files Hands_off_control_HLS/l1_control_sol1.cpp
add_files Hands_off_control_HLS/ltris.h
add_files -tb Hands_off_control_HLS/hoc_sim_tb.cpp -cflags "-Wno-unknown-pragmas"
open_solution "solution1" -flow_target vivado
set_part {xc7z020-clg400-3}
create_clock -period 10 -name default
config_export -format ip_catalog -output /home/robsta/Desktop/HLS_Projects/Hoc_final_results -rtl verilog
source "./Hands_off_control_HLS/solution1/directives.tcl"
csim_design
csynth_design
cosim_design
export_design -rtl verilog -format ip_catalog -output /home/robsta/Desktop/HLS_Projects/Hoc_final_results
