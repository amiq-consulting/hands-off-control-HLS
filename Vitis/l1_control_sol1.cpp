#include "l1_control.h"

void hoc_module(hls::stream<io_variable_axis_wrapper> &in_state, hls::stream<io_variable_axis_wrapper> &out_port) {

	#pragma HLS INTERFACE mode=ap_ctrl_none port=return

	#pragma HLS INTERFACE mode=axis register_mode=both port=in_state register
	#pragma HLS INTERFACE mode=axis register_mode=both port=out_port register

	ADMM_variable colo_min_time;
	ap_int<32> temp_min_time = 0;

	// ADMM Variables

	ADMM_variable reg_aux_var_1 = (ADMM_variable)0;
	ADMM_variable reg_aux_var_2 = (ADMM_variable)0;

	ADMM_variable reg_vec_b[max_size] = {(ADMM_variable)0};

	ADMM_variable reg_vec_x[2 * max_size] = {(ADMM_variable)0};	 					// Array that stores control signals in first half and slack variables in second half
	ADMM_variable reg_vec_y[3 * max_size + 2] = {(ADMM_variable)0};					// Array that stores auxiliary ADMM variable 1
	ADMM_variable reg_vec_z[3 * max_size + 2] = {(ADMM_variable)0};					// Array that stores auxiliary ADMM variable 2

	ADMM_variable reg_vec_x_next[2 * max_size] = {(ADMM_variable)0};				// Array that stores shifted x vector
	ADMM_variable reg_vec_y_next[3 * max_size + 2] = {(ADMM_variable)0};			// Array that stores shifted y vector
	ADMM_variable reg_vec_z_next[3 * max_size + 2] = {(ADMM_variable)0};			// Array that stores shifted z vector

	ADMM_variable reg_vec_b1[max_size] = {(ADMM_variable)0};						// Array that stores right hand side of KKT system - first n_L1
	ADMM_variable reg_vec_b2[max_size] = {(ADMM_variable)0};						// Array that stores right hand side of KKT system - second n_L1
	ADMM_variable reg_vec_b3[max_size] = {(ADMM_variable)0};						// Array that stores right hand side of KKT system - third n_L1
	ADMM_variable reg_vec_b4[max_size] = {(ADMM_variable)0};						// Array that stores right hand side of KKT system - fourth n_L1
	ADMM_variable reg_vec_b5[max_size] = {(ADMM_variable)0};						// Array that stores right hand side of KKT system - fifth n_L1
	ADMM_variable reg_vec_b6[2] = {(ADMM_variable)0};								// Array that stores right hand side of KKT system - last two elements

	ADMM_variable reg_vec_u[max_size] = {(ADMM_variable)0};							// Array that stores control signals
	ADMM_variable reg_vec_u_temp[max_size] = {(ADMM_variable)0};
	ADMM_variable reg_vec_t[max_size] = {(ADMM_variable)0};							// Array that stores slack variables

	ADMM_variable reg_vec_v_k_1[3 * max_size + 2] = {(ADMM_variable)0};

	ADMM_variable reg_vec_temp_z_next[3 * max_size + 2] = {(ADMM_variable)0};

	ADMM_variable reg_alpha = (ADMM_variable)alfa;
	ADMM_variable reg_sigma = (ADMM_variable)sigma;
	ADMM_variable reg_inv_rho = (ADMM_variable)inv_rho;
	ADMM_variable reg_rho = (ADMM_variable)rho;
	ADMM_variable reg_inv_sigma_3_rho = (ADMM_variable)(1.0 / (sigma + 3 * rho));
	ADMM_variable reg_inv_r = (ADMM_variable)(1.0 / r);

	ADMM_variable reg_mat_E_line_1[max_size] = {(ADMM_variable)0};
	ADMM_variable reg_mat_E_line_2[max_size] = {(ADMM_variable)0};

//	ADMM_variable reg_mat_L_tr[max_size][max_size] = {(ADMM_variable)0};
	ADMM_variable reg_mat_L_inf_tr[(max_size + 1) * (max_size / 2)] = {(ADMM_variable)0};

	//

	state_variable reg_Ka = (state_variable)Ka;
	state_variable reg_Kb = (state_variable)Kb;

	state_variable reg_Ka_d = (state_variable)Ka_d;
	state_variable reg_Kb1_d = (state_variable)Kb1_d;
	state_variable reg_Kb2_d = (state_variable)Kb2_d;
	state_variable reg_Ka_d_Kb1_d = (state_variable)(Ka_d * Kb1_d);

	state_variable reg_inv_Ka_Kb = (state_variable)(1.0 / (Ka * Kb));
	state_variable reg_inv_Kb = (state_variable)(1.0 / Kb);

	state_variable reg_inv_Ts = (state_variable)(1.0 / Ts);

	io_variable_axis_wrapper io_temp_in;
	io_variable_axis_wrapper io_temp_out;

	io_variable io_q_state;
	io_variable io_w_state;

	state_variable q_state = (state_variable)0;
	state_variable w_state = (state_variable)0;

	state_variable abs_w_state = (state_variable)0; 			// Store absolute value of w0
	state_variable aux_compare_variable = (state_variable)0;	// Store value required for comparison
	state_variable temp_w0_kb_division = (state_variable)0;		// Store value of w0/kb
	state_variable temp_w0_kb_squared = (state_variable)0;		// Store value of (w0/kb)^2
	state_variable temp_4_q0_inv_ka_kb = (state_variable)0; 	// Store value of 4*q0/(ka*kb)
	state_variable temp_sqrt = (state_variable)0;				// Store value to be square rooted

	state_variable local_min_time = (state_variable)0;
	state_variable reg_temp = (state_variable)0;

	in_state.read(io_temp_in);
	io_q_state.range() = io_temp_in.data;
	q_state = io_q_state;

	in_state.read(io_temp_in);
	io_w_state.range() = io_temp_in.data;
	w_state = io_w_state;

	// Computing minimum time

	abs_w_state = hls::fabs(w_state); // Absolute value of w0

	temp_w0_kb_division = w_state * reg_inv_Kb;															// Save constant value 1
	temp_w0_kb_squared = temp_w0_kb_division * temp_w0_kb_division;										// Save constant value 2

	aux_compare_variable = ((state_variable)-0.5) * abs_w_state * temp_w0_kb_division * reg_Ka;			// Save value to compare with

	temp_4_q0_inv_ka_kb = 4 * q_state * reg_inv_Ka_Kb;													//


	if (q_state > aux_compare_variable) {
		temp_sqrt = temp_4_q0_inv_ka_kb;
		local_min_time = temp_w0_kb_division;
	} else {
		temp_sqrt = -temp_4_q0_inv_ka_kb;
		local_min_time = -temp_w0_kb_division;
	}

	temp_sqrt = temp_sqrt + 2 * temp_w0_kb_squared;
	local_min_time = local_min_time + hls::sqrt(temp_sqrt);

	reg_temp = (ADMM_variable)reg_inv_r * (ADMM_variable)local_min_time * (ADMM_variable)reg_inv_Ts;						// Transform into samples

	colo_min_time = hls::ceil(reg_temp);

	if (colo_min_time <= (ADMM_variable)20)
		colo_min_time = 21;
	temp_min_time = colo_min_time;

	// Solving L1 Problem
	if ((colo_min_time) < max_size && (colo_min_time > (ADMM_variable)(20))) {

		chol_decomp<ADMM_variable, max_size>(reg_mat_L_inf_tr, temp_min_time);

		reg_aux_var_1 = 0;
//		compute_mat_E_loop: for (int i = 0; i < colo_min_time; i++) {
//			reg_mat_E_line_1[i] = reg_Kb1_d;
//			reg_mat_E_line_2[i] = (colo_min_time - 1 - reg_aux_var_1) * reg_Ka_d * reg_Kb1_d + reg_Kb2_d;
//			reg_aux_var_1 += 1;
//		}
		compute_mat_E_loop: for (int i = 0; i < max_size; i++) {
			if (i < colo_min_time) {
				reg_mat_E_line_1[i] = reg_Kb1_d;
				reg_mat_E_line_2[i] = (colo_min_time - 1 - reg_aux_var_1) * reg_Ka_d_Kb1_d + reg_Kb2_d;
			}
			reg_aux_var_1 += 1;
		}


		// Solving L1 Problem
		ADMM_loop: for (int k = 0; k < n_of_iterations; k++) {
			// Assigning bi vectors
			b_vector_assigning_loop: for (int i = 0; i < max_size; i++) {
				// b = [sigma * x_k - q; z_k - y_k / rho]
				// b1 = first n_L1 elements of b
				// b2 = second n_L1 elements of b
				// b3 = third n_L1 elements of b
				// b4 = fourth n_L1 elements of b
				// b5 = fifth n_L1 elements of b
				// b6 = last two elements of b

				reg_vec_b1[i] = reg_sigma * reg_vec_x[i];
				reg_vec_b2[i] = reg_sigma * reg_vec_x[max_size + i] - (ADMM_variable)1;

				reg_vec_b3[i] = reg_vec_z[i] - reg_inv_rho * reg_vec_y[i];
				reg_vec_b4[i] = reg_vec_z[i + max_size] - reg_inv_rho * reg_vec_y[i + max_size];
				reg_vec_b5[i] = reg_vec_z[i + max_size + max_size] - reg_inv_rho * reg_vec_y[i + max_size + max_size];
			}

			reg_vec_b6[0] = reg_vec_z[3 * max_size] - reg_inv_rho * reg_vec_y[3 * max_size];
			reg_vec_b6[1] = reg_vec_z[3 * max_size + 1] - reg_inv_rho * reg_vec_y[3 * max_size + 1];

			// Computing vector to be used for linear solver
			u_t_vector_computation_loop: for (int i = 0; i < max_size; i++) {
				reg_vec_u[i] = reg_mat_E_line_1[i] * reg_vec_b6[0] + reg_mat_E_line_2[i] * reg_vec_b6[1] +
						reg_vec_b3[i] + reg_vec_b4[i] + reg_vec_b1[i] * reg_inv_rho;
				reg_vec_t[i] = (reg_rho * (-reg_vec_b3[i] + reg_vec_b4[i] + reg_vec_b5[i]) + reg_vec_b2[i])
										* reg_inv_sigma_3_rho;
			}

			ltris<ADMM_variable, max_size>(reg_mat_L_inf_tr, reg_vec_u, reg_vec_u_temp, temp_min_time);

			u_u_temp_vector_transfer_loop_1: for (int i = 0; i < max_size; i++) {
				reg_vec_u[i] = reg_vec_u_temp[i];
			}

			utris<ADMM_variable, max_size>(reg_mat_L_inf_tr, reg_vec_u, reg_vec_u_temp, temp_min_time);

			u_u_temp_vector_transfer_loop_2: for (int i = 0; i < max_size; i++) {
				reg_vec_u[i] = reg_vec_u_temp[i];
			}

			reg_aux_var_1 = (ADMM_variable)0;
			reg_aux_var_2 = (ADMM_variable)0;

			// Computing vector v
			v_vector_computation_loop: for (int i = 0; i < max_size; i++) {
				reg_vec_v_k_1[i] 				= (reg_vec_u[i] - reg_vec_t[i] - reg_vec_b3[i]) * reg_rho;
				reg_vec_v_k_1[i + max_size] 	= (reg_vec_u[i] + reg_vec_t[i] - reg_vec_b4[i]) * reg_rho;
				reg_vec_v_k_1[i + 2 * max_size] = (reg_vec_t[i] - reg_vec_b5[i]) * reg_rho;
			}

//			aux_E_b6_multiplication_loop: for (int i = 0; i < colo_min_time; i++) {
//				reg_aux_var_1 += reg_vec_u[i];
//				reg_aux_var_2 += (colo_min_time - i - 1) * reg_vec_u[i];
//			}
			aux_E_b6_multiplication_loop: for (int i = 0; i < max_size; i++) {
				if (i < colo_min_time) {
					reg_aux_var_1 += reg_vec_u[i];
					reg_aux_var_2 += (colo_min_time - i - 1) * reg_vec_u[i];
				}
			}

			reg_vec_v_k_1[3 * max_size] = (reg_Kb1_d * reg_aux_var_1 - reg_vec_b6[0]) * reg_rho;
			reg_vec_v_k_1[3 * max_size + 1] = (reg_Ka_d_Kb1_d * reg_aux_var_2 + reg_Kb2_d * reg_aux_var_1 - reg_vec_b6[1]) * reg_rho;

			// Computing temp_z_k_1 vector
			z_next_computation_loop: for (int i = 0; i < 3 * max_size + 2; i++) {
				reg_vec_temp_z_next[i] = reg_vec_z[i] + (reg_vec_v_k_1[i] - reg_vec_y[i]) * reg_inv_rho;
			}

			// Updating x_k vector
			x_next_computation_loop: for (int i = 0; i < max_size; i++) {
				reg_vec_x_next[i] 				= reg_alpha * reg_vec_u[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_x[i];
				reg_vec_x_next[i + max_size] 	= reg_alpha * reg_vec_t[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_x[i + max_size];
			}

//			// Updating z_k vector
//			z_update_loop1: for (int i = 0; i < max_size; i++) {
//				reg_aux_var_1 = reg_alpha * reg_vec_temp_z_next[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_z[i] + reg_vec_y[i] * reg_inv_rho;
//				if (reg_aux_var_1 > 0)
//					reg_aux_var_1 = 0;
//				else if (reg_aux_var_1 < -2)
//					reg_aux_var_1 = -2;
//				reg_vec_z_next[i] = reg_aux_var_1;
//			}
//
//			z_update_loop2: for (int i = max_size; i < 2 * max_size; i++) {
//				reg_aux_var_1 = reg_alpha * reg_vec_temp_z_next[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_z[i] + reg_vec_y[i] * reg_inv_rho;
//				if (reg_aux_var_1 > 2)
//					reg_aux_var_1 = 2;
//				else if (reg_aux_var_1 < 0)
//					reg_aux_var_1 = 0;
//				reg_vec_z_next[i] = reg_aux_var_1;
//			}
//
//			z_update_loop3: for (int i = 2 * max_size; i < 3 * max_size; i++) {
//				reg_aux_var_1 = reg_alpha * reg_vec_temp_z_next[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_z[i] + reg_vec_y[i] * reg_inv_rho;
//				if (reg_aux_var_1 > 1)
//					reg_aux_var_1 = 1;
//				else if (reg_aux_var_1 < 0)
//					reg_aux_var_1 = 0;
//				reg_vec_z_next[i] = reg_aux_var_1;
//			}

			z_update_loop : for(int i = 0; i < 3 * max_size; i++) {
				reg_vec_z_next[i] = reg_alpha * reg_vec_temp_z_next[i] + ((ADMM_variable)1 - reg_alpha) * reg_vec_z[i] + reg_vec_y[i] * reg_inv_rho;
			}

			z_project_loop : for(int i = 0; i < max_size; i++) {
				if (reg_vec_z_next[i] > 0)
					reg_vec_z_next[i] = 0;
				else if (reg_vec_z_next[i] < -2)
					reg_vec_z_next[i] = -2;

				if (reg_vec_z_next[i + max_size] > 2)
					reg_vec_z_next[i + max_size] = 2;
				else if (reg_vec_z_next[i + max_size] < 0)
					reg_vec_z_next[i + max_size] = 0;

				if (reg_vec_z_next[i + 2 * max_size] > 1)
					reg_vec_z_next[i + 2 * max_size] = 1;
				else if (reg_vec_z_next[i + 2 * max_size] < 0)
					reg_vec_z_next[i + 2 * max_size] = 0;
			}

			reg_vec_z_next[3 * max_size] = -(ADMM_variable)w_state;
			reg_vec_z_next[3 * max_size + 1] =  -colo_min_time * reg_Ka_d * (ADMM_variable)w_state - (ADMM_variable)q_state;

			// Updating y vector
			y_update_loop: for (int i = 0; i < 3 * max_size + 2; i++) {
				reg_vec_y_next[i] = reg_vec_y[i] + reg_rho * (reg_alpha * reg_vec_temp_z_next[i] + (1 - reg_alpha) * reg_vec_z[i] - reg_vec_z_next[i]);
			}

			// Shifting vector x
			x_shift_loop: for (int i = 0; i < 2 * max_size; i++) {
				reg_vec_x[i] = reg_vec_x_next[i];
			}

			// Shifting vector y and z
			y_z_shift_loop: for (int i = 0; i < 3 * max_size + 2; i++){
				reg_vec_y[i] = reg_vec_y_next[i];
				reg_vec_z[i] = reg_vec_z_next[i];
			}
		}
	}

	io_variable out_aux_var;

	out_aux_var = local_min_time;
	io_temp_out.data = out_aux_var.range();
	out_port.write(io_temp_out);

	out_aux_var = colo_min_time;
	io_temp_out.data = out_aux_var.range();
	out_port.write(io_temp_out);

	for (int i = 0; i < max_size; i++) {
		out_aux_var = reg_vec_x[i];
		io_temp_out.data = out_aux_var.range();
		out_port.write(io_temp_out);
	}
}
