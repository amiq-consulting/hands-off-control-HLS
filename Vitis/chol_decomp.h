#ifndef CHOL_DECOMP_H
#define CHOL_DECOMP_H

template <typename T, int N>
void chol_decomp(T L[(N / 2) * (N + 1)], int L1_problem_size) {
#pragma HLS INLINE

	static T reg_Kb_d_sq_sum 		= (T)(Kb1_d * Kb1_d + Kb2_d * Kb2_d);
	static T reg_Ka_d_sq_Kb1_d_sq 	= (T)(Ka_d * Ka_d * Kb1_d * Kb1_d);
	static T reg_Ka_d_Kb1_d_Kb2_d 	= (T)(Ka_d * Kb1_d * Kb2_d);
	static T reg_sigma_inv_rho 		= (T)(sigma / rho);

	int reg_int_aux_index_1 = 0;
	int reg_int_aux_index_2 = 0;

	T reg_L1_problem_size = (T)L1_problem_size;

	T reg_alfa = (T)0;
	T reg_inv_alfa = (T)0;
	T reg_aux = (T)0;
	T reg_temp = (T)0;
	T reg_constant_diag_value = (T)constant_diag_value;

	T reg_store_1 = (T)0;
	T reg_store_2 = (T)0;

	T reg_vec_store[max_size] = {(T)0};
	T reg_vec_rez[max_size] = {(T)0};

	reg_alfa = (L1_problem_size - 1) * (L1_problem_size - 1) * reg_Ka_d_sq_Kb1_d_sq +
			(L1_problem_size + L1_problem_size - 2) * reg_Ka_d_Kb1_d_Kb2_d + reg_Kb_d_sq_sum + (T)(2 + reg_sigma_inv_rho);
	reg_alfa = hls::sqrt(reg_alfa);
	L[0] = reg_alfa;
	reg_inv_alfa = (T)1.0 / reg_alfa;

	reg_int_aux_index_1 = 1;
	chol_first_col_computation_loop: for (int i = 1; i < max_size; i++) {
//#pragma HLS PIPELINE
		reg_aux = 0;
		if (i < L1_problem_size) {
			reg_aux = ((L1_problem_size - ((T)i + 1)) * (L1_problem_size - 1) * reg_Ka_d_sq_Kb1_d_sq +
					(L1_problem_size - ((T)i + 1) + L1_problem_size - 1) * reg_Ka_d_Kb1_d_Kb2_d + reg_Kb_d_sq_sum);
		}
			L[reg_int_aux_index_1] = reg_aux * reg_inv_alfa;;

			reg_int_aux_index_1 += i + 1;
	}

	chol_col_iteration_loop: for (int k = 1; k < max_size; k++) {
		if (k < L1_problem_size) {
			reg_alfa = (L1_problem_size - (T)k - 1) * (L1_problem_size - (T)k - 1) * reg_Ka_d_sq_Kb1_d_sq +
						(L1_problem_size - (T)k + L1_problem_size - (T)k - 2) * reg_Ka_d_Kb1_d_Kb2_d + reg_Kb_d_sq_sum + (T)(2 + reg_sigma_inv_rho);

			reg_int_aux_index_1 = (k * (k + 1)) / 2;
			chol_diag_elem_computation_loop: for (int s = 0; s < max_size; s++) {
#pragma HLS PIPELINE
				if (s < k)
//					reg_alfa -= L[(k * (k + 1)) / 2 + s] * L[(k * (k + 1)) / 2 + s];
					reg_alfa -= L[reg_int_aux_index_1 + s] * L[reg_int_aux_index_1 + s];
			}

			reg_alfa = hls::sqrt(reg_alfa);
			reg_inv_alfa = (T)1.0 / reg_alfa;
//			L[(k * (k + 1)) / 2 + k] = reg_alfa;
			L[reg_int_aux_index_1 + k] = reg_alfa;

			chol_col_store_loop: for (int s = 0; s < max_size; s++) {
#pragma HLS PIPELINE
				if (s < k)
//					reg_vec_store[s] = L[(k * (k + 1)) / 2 + s];
					reg_vec_store[s] = L[reg_int_aux_index_1 + s];
			}

			chol_col_computation_loop: for (int i = 0; i < max_size; i++) {
				if ((i + k + 1) < L1_problem_size) {
					reg_aux = ((L1_problem_size - (T)i - 2 - (T)k) * (L1_problem_size - (T)k - 1) * reg_Ka_d_sq_Kb1_d_sq +
							(L1_problem_size - (T)i + L1_problem_size - (T)k - 2 - (T)k - 1) * reg_Ka_d_Kb1_d_Kb2_d + reg_Kb_d_sq_sum);

					reg_int_aux_index_2 = ((i + k + 1) * (i + k + 2)) / 2;
					chol_col_elem_computation_loop: for (int s = 0; s < max_size; s++) {
#pragma HLS PIPELINE
						if (s < k) {
//							reg_aux -= L[((i + k + 1) * (i + k + 2)) / 2 + s] * reg_vec_store[s];
							reg_aux -= L[reg_int_aux_index_2 + s] * reg_vec_store[s];
						}
					}

//					reg_vec_rez[i + k + 1] = reg_aux / reg_alfa;
					reg_vec_rez[i + k + 1] = reg_aux * reg_inv_alfa;
				}
			}

			reg_int_aux_index_2 = ((k + 1) * (k + 2)) / 2 + k;
			chol_col_results_store_loop: for (int i = 0; i < max_size; i++) {
#pragma HLS PIPELINE
				if ((i + k + 1) < L1_problem_size)
//					L[((i + k + 1) * (i + k + 2)) / 2 + k] = reg_vec_rez[i + k + 1];
					L[reg_int_aux_index_2] = reg_vec_rez[i + k + 1];

				reg_int_aux_index_2 += i + k + 2;
			}
		}
	}

//	chol_rest_of_matrix_computation_loop: for (int k = L1_problem_size; k < max_size; k++) {
//		L[(k * (k + 1)) / 2 + k] = reg_constant_diag_value;
//	}

	reg_int_aux_index_1 = 0;
	chol_rest_of_matrix_computation_loop: for (int qq = 0; qq < max_size; qq++) {
		if (qq >= reg_L1_problem_size)
			L[reg_int_aux_index_1] = reg_constant_diag_value;

		reg_int_aux_index_1 += qq + 2;
	}
}

#endif
