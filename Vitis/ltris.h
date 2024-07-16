#ifndef LTRIS_H
#define LTRIS_H

template <typename T, int N>
void ltris(T L[(N * (N + 1)) / 2], T b[N], T x[N], int L1_problem_size) {
#pragma HLS INLINE

	T reg_aux_var = 0;
	ap_uint<24> reg_aux_index = 0;

	ltris_front_substitution_loop: for (int i = 0; i < max_size; i++) {
//#pragma HLS PIPELINE II=401
		reg_aux_var = b[i];

		reg_aux_index = (i * (i + 1)) / 2;
		if (i < L1_problem_size) {
			ltris_current_element_computation_loop: for (int j = 0; j < max_size; j++) {
				if (j < i)
					reg_aux_var -= L[reg_aux_index + j] * x[j];
			}
		}

//		double aux = L[reg_aux_index + i];
//		if (i > 250)
//			printf("Puna siua %d %f\n", reg_aux_index + i, aux);

		x[i] = reg_aux_var / L[reg_aux_index + i];
	}
//		reg_aux_index = (i * (i + 1)) / 2;
//		ltris_current_element_computation_loop: for (int j = 0; j < max_size; j++) {
////#pragma HLS PIPELINE
//			if ((j < i) && (i < L1_problem_size))
//				reg_aux_var -= L[reg_aux_index + j] * x[j];
//		}
//


//
//	ltris_front_substitution_loop: for (int i = 0; i < L1_problem_size; i++) {
//		reg_aux_var = b[i];
//
//		ltris_current_element_computation_loop: for (int j = 0; j < i; j++) {
//			reg_aux_var -= L[(i * (i + 1)) / 2 + j] * x[j];
//		}
//
//		x[i] = reg_aux_var / L[(i * (i + 1)) / 2 + i];
//	}
//
//	ltris_rest_of_vector_loop: for (int i = L1_problem_size; i < N; i++) {
//		x[i] /= L[(i * (i + 1)) / 2 + i];
//	}
}

template <typename T, int N>
void utris(T U[(N * (N + 1)) / 2], T b[N], T x[N], int L1_problem_size) {
#pragma HLS INLINE

	T reg_aux_var = 0;
	ap_uint<24> reg_aux_index = 0;

	utris_back_substitution_loop: for (int i = max_size - 1; i >= 0; i--) {
		reg_aux_var = b[i];

//		reg_aux_index = ((i + 1) * (i + 2)) / 2 + i;////////////////////////////////////
//		if (i < L1_problem_size) {
//			utris_current_element_computation_loop : for (int j = 0; j < max_size; j++) {
//				if (j < (L1_problem_size - i - 1))
//					reg_aux_var -= U[reg_aux_index] * x[i + j + 1];
//
//				reg_aux_index += i + j + 2;
//			}
//		}

//		reg_aux_index = ((i + 1) * (i + 2)) / 2 + i;////////////////////////////////////
		reg_aux_index = 0;
		if (i < L1_problem_size) {
			utris_current_element_computation_loop : for (int j = 0; j < max_size; j++) {
				if ((j >= (i + 1)) && (j < L1_problem_size))
					reg_aux_var -= U[reg_aux_index + i] * x[j];

//				reg_aux_index += i + j + 2;
				reg_aux_index += j + 1;
			}
		}

		x[i] = reg_aux_var / U[(i * (i + 1)) / 2 + i];
	}

//	utris_back_substitution_loop: for (int i = max_size - 1; i >= 0; i--) {
////#pragma HLS PIPELINE II=401
//		reg_aux_var = b[i];
//
//		reg_aux_index = ((i + 1) * (i + 2)) / 2 + i;
//		utris_current_element_computation_loop : for (int j = 0; j < max_size; j++) {
////#pragma HLS PIPELINE
//			if ((i < L1_problem_size) && ((j + i + 1) < L1_problem_size))
//				reg_aux_var -= U[reg_aux_index] * x[i + j + 1];
//
//			reg_aux_index += i + j + 2;
//		}
//
//		x[i] = reg_aux_var / U[(i * (i + 1)) / 2 + i];
//	}

//	utris_rest_of_vector_loop: for (int i = N - 1; i >= L1_problem_size; i--) {
//		x[i] /= U[(i * (i + 1)) / 2 + i];
//	}
//
//	utris_back_substitution_loop: for (int i = L1_problem_size - 1; i >= 0; i--) {
//		reg_aux_var = b[i];
//
//		utris_current_element_computation_loop: for (int j = i + 1; j < L1_problem_size; j++) {
//			reg_aux_var -= U[(j * (j + 1)) / 2 + i] * x[j];
//		}
//
//		x[i] = reg_aux_var / U[(i * (i + 1)) / 2 + i];
//	}

}

#endif
