open_project Hands_off_control_HLS
set_top hoc_module
add_files chol_decomp.h
add_files l1_control.h
add_files l1_control_sol1.cpp
add_files ltris.h
add_files -tb hoc_sim_tb.cpp -cflags "-Wno-unknown-pragmas"
open_solution "solution1" -flow_target vivado
set_part {xc7z020-clg400-3}
create_clock -period 10 -name default
csim_design
csynth_design
cosim_design
