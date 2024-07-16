#include "l1_control.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <random>

void min_time_sw(double &state_q0, double &state_w0, double &min_time);
void l1_control_sw(double &state_q0, double &state_w0, int L1_horizon, double *u_opt);

double mmin(double a, double b) {
	return (a < b)? a : b;
}

double mmax(double a, double b) {
	return (a > b)? a : b;
}

int main() {
	int n_of_steps = 30;

	int global_it = 0;

	FILE *state_file;
	FILE *timing_file;
	FILE *software_output_file;
	FILE *hardware_output_file;

	srand((unsigned int)time(NULL));

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());

    std::normal_distribution<double> dist(0, 0.003);

//    int i;
//    float sample;
//    for(i = 0; i < 1000; ++i)
//    {
//        // instance of class std::normal_distribution with specific mean and stddev
//        std::normal_distribution<float> d(mean[i], stddev[i]);
//
//        // get random number with normal distribution using gen as random source
//        sample = d(gen);
//
//        // profit
//        do_something_with_this_value(sample);
//    }

	io_variable aux_var;
	io_variable_axis_wrapper local_state;
	hls::stream<io_variable_axis_wrapper> in_state;
	hls::stream<io_variable_axis_wrapper> min_time_hw;

	double state_q0 = 3;
	double state_w0 = 0.5;

	double s_q0 = state_q0;
	double s_w0 = state_w0;

	double state_x_curr = state_q0;
	double state_x_dot_curr = state_w0;

	double sw_min_time[n_of_steps];
	double hw_min_time[n_of_steps];

	double hw_u_opt[10000];
	double sw_u_opt[10000];

	double state_x_vec[10000];
	double state_x_dot_vec[10000];

	double hw_temp_u[max_size];

	double l1_horizon = 0;


	for (int i = 0; i < n_of_steps; i++) {
		// Calling the software function for minimum time
		min_time_sw(s_q0, s_w0, sw_min_time[i]);

		// Calling the hardware function
		aux_var = (io_variable)s_q0;
		local_state.data = aux_var.range();
		in_state.write(local_state);

		aux_var = (io_variable)s_w0;
		local_state.data = aux_var.range();
		in_state.write(local_state);

		hoc_module(in_state, min_time_hw);

		// Read the minimum time
		min_time_hw.read(local_state);
		aux_var.range() = local_state.data;
		hw_min_time[i] = aux_var;

		// Read the L1 horizon
		min_time_hw.read(local_state);
		aux_var.range() = local_state.data;
		l1_horizon = aux_var;

		printf("Software L1 horizon: %f, Hardware L1 horizon: %f \n", ceil(sw_min_time[i] / (Ts * r)), l1_horizon);

		// Check if horizon does not exceed maximum length and is also not too short
		if ((l1_horizon >= max_size) || (l1_horizon <= 20)) {

			for (int j = 0; j < max_size; j++)
				min_time_hw.read(local_state);

			state_x_vec[global_it] = state_x_curr;
			state_x_dot_vec[global_it] = state_x_dot_curr;
			hw_u_opt[global_it] = 0;
			sw_u_opt[global_it] = 0;

			global_it++;

			state_x_curr += Ka_d * state_x_dot_curr;
//			state_x_dot_curr += 0.005 * (((double)rand() / (float)(RAND_MAX)) * 2 - 1);
			state_x_dot_curr += mmax(-0.01, mmin(0.01, dist(gen)));

			printf("BUNA l1 horizon e %f\n", l1_horizon);

			continue;
		}

		// Receive hardware optimal control
		for (int j = 0; j < max_size; j++) {
			min_time_hw.read(local_state);
			aux_var.range() = local_state.data;
			hw_temp_u[j] = aux_var;
			hw_temp_u[j] = round(hw_temp_u[j]);
		}

		// Save software optimal control
		l1_control_sw(s_q0, s_w0, ceil(sw_min_time[i] / (Ts * r)), sw_u_opt + global_it);

		// Apply hardware optimal contrl
		for (int j = 0; j < l1_horizon; j++) {
			state_x_vec[global_it] = state_x_curr;
			state_x_dot_vec[global_it] = state_x_dot_curr;

			state_x_curr += Ka_d * state_x_dot_curr + Kb2_d * hw_temp_u[j];
//			state_x_dot_curr += Kb1_d * hw_temp_u[j] + 0.008 * (((double)rand() / (float)(RAND_MAX)) * 2 - 1);
			state_x_dot_curr += Kb1_d * hw_temp_u[j] + mmax(-0.01, mmin(0.01, dist(gen)));

			hw_u_opt[global_it++] = hw_temp_u[j];
		}

		s_q0 = state_x_curr;
		s_w0 = state_x_dot_curr;
	}

	// Saving output data to separate files

	// State variables
	state_file = fopen("state_file.dat", "w");

	for (int i = 0; i < global_it; i++) {
		fprintf(state_file, "%.5f\t", state_x_vec[i]);
		fprintf(state_file, "%.5f\n", state_x_dot_vec[i]);
	}

	fclose(state_file);

//	// Timings
//	timing_file = fopen("timing_file.dat", "w");
//
//	fprintf(timing_file, "Hardware Minimum TIme = %f seconds\n", hw_min_time[0]);
//	fprintf(timing_file, "Software Minimum Time = %f seconds\n", sw_min_time[0]);
//
//	fprintf(timing_file, "Software L1 horizon: %d\n Hardware L1 horizon: %f \n", L1_horizon, aux);
//
//	fclose(timing_file);

	// Software output file
	software_output_file = fopen("software_output_file.dat", "w");

	for (int i = 0; i < global_it; i++) {
		fprintf(software_output_file, "%.5f\n", sw_u_opt[i]);
	}

	fclose(software_output_file);

	// Hardware output file
	hardware_output_file = fopen("hardware_output_file.dat", "w");

	for (int i = 0; i < global_it; i++) {
		fprintf(hardware_output_file, "%.5f\n", hw_u_opt[i]);
	}

	fclose(hardware_output_file);

	// Done outputting to files

	return 0;
}

void min_time_sw(double &state_q0, double &state_w0, double &min_time) {
	double aux_variable = -0.5 * state_w0 * fabs(state_w0) * Ka / Kb;

	if (state_q0 > aux_variable)
		min_time = state_w0 / Kb + sqrt(2 * pow(state_w0 / Kb, 2) + 4 * state_q0 / (Ka * Kb));
	else
		min_time = -state_w0 / Kb + sqrt(2 * pow(state_w0 / Kb, 2) - 4 * state_q0 / (Ka * Kb));

}

void l1_control_sw(double &state_q0, double &state_w0, int L1_horizon, double *u_opt) {
	double x_k[2 * L1_horizon] = {0};
	double x_k_1[2 * L1_horizon] = {0};

	double z_k[3 * L1_horizon + 2] = {0};
	double z_k_1[3 * L1_horizon + 2] = {0};

	double y_k[3 * L1_horizon + 2] = {0};
	double y_k_1[3 * L1_horizon + 2] = {0};

	double b1[L1_horizon] = {0};
	double b2[L1_horizon] = {0};
	double b3[L1_horizon] = {0};
	double b4[L1_horizon] = {0};
	double b5[L1_horizon] = {0};
	double b6[2] = {0};

	double u[max_size] = {0};
	double u_temp[max_size] = {0};
	double t[L1_horizon] = {0};

	double v1[L1_horizon] = {0};
	double v2[L1_horizon] = {0};
	double v3[L1_horizon] = {0};
	double v4[2] = {0};
	double v_k_1[3 * L1_horizon + 2] = {0};

	double kkt_mat[L1_horizon * L1_horizon] = {0};
	double L_mat[max_size][max_size] = {0};
	double L_mat_inf_tr[max_size * (max_size + 1) / 2] = {0};

	double aux_acc_1, aux_acc_2 = {0};

	double temp_z_k_1[3 * L1_horizon + 2] = {0};

	int info;

	for (int row = 0; row < L1_horizon; row++) {
		for (int col = 0; col < L1_horizon; col++) {
			kkt_mat[row * L1_horizon + col] = Kb1_d * Kb1_d + Kb2_d * Kb2_d + (L1_horizon - row - 1) * (L1_horizon - col - 1) * Ka_d * Ka_d * Kb1_d * Kb1_d +
					(2 * L1_horizon - row - col - 2) * Ka_d * Kb1_d * Kb2_d;
//			kkt_mat[col * L1_horizon + row] = kkt_mat[row * L1_horizon + col];

			if (row == col)
				kkt_mat[row * L1_horizon + col] += 2 + sigma / rho;
		}
	}

	chol_decomp<double, max_size>(L_mat_inf_tr, L1_horizon);

//	for (int i = 0; i < max_size; i++)
//		printf("Elem %d la pozitia %d: %f\n", i, (i * (i + 1)) / 2 + i, L_mat_inf_tr[(i * (i + 1)) / 2 + i]);

	for (int i = 0; i < n_of_iterations; i++) {
		for (int j = 0; j < L1_horizon; j++) {
			b1[j] = sigma * x_k[j];
			b2[j] = sigma * x_k[j + L1_horizon] - 1;

			b3[j] = z_k[j] - y_k[j] / rho;
			b4[j] = z_k[j + L1_horizon] - y_k[j + L1_horizon] / rho;
			b5[j] = z_k[j + 2 * L1_horizon] - y_k[j + 2 * L1_horizon] / rho;
		}

		b6[0] = z_k[3 * L1_horizon] - y_k[3 * L1_horizon] / rho;
		b6[1] = z_k[3 * L1_horizon + 1] - y_k[3 * L1_horizon + 1] / rho;

		for (int j = 0; j < L1_horizon; j++) {
			u[j] = Kb1_d * b6[0] + ((L1_horizon - j - 1) * Ka_d * Kb1_d + Kb2_d) * b6[1] + b3[j] + b4[j] + b1[j] / rho;
			t[j] = (rho * (-b3[j] + b4[j] + b5[j]) + b2[j]) / (sigma + 3 * rho);
		}

//		xf_orig::solver_orig::polinearsolver<double, max_size, 10>(L1_horizon, kkt_mat, 1, u, L1_horizon, 1, info);

		ltris<double, max_size>(L_mat_inf_tr, u, u_temp, L1_horizon);

		for (int j = 0; j < max_size; j++)
			u[j] = u_temp[j];

		utris<double, max_size>(L_mat_inf_tr, u, u_temp, L1_horizon);

		for (int j = 0; j < max_size; j++)
			u[j] = u_temp[j];

		aux_acc_1 = 0;
		aux_acc_2 = 0;

		for (int j = 0; j < L1_horizon; j++) {
			v_k_1[j] 					= rho * (u[j] - t[j] - b3[j]);
			v_k_1[j + L1_horizon] 		= rho * (u[j] + t[j] - b4[j]);
			v_k_1[j + 2 * L1_horizon] 	= rho * (t[j] - b5[j]);

			aux_acc_1 += u[j];
			aux_acc_2 += u[j] * (L1_horizon - j - 1);
		}

		v_k_1[3 * L1_horizon] 			= rho * (Kb1_d * aux_acc_1 - b6[0]);
		v_k_1[3 * L1_horizon + 1] 		= rho * (Ka_d * Kb1_d * aux_acc_2 + Kb2_d * aux_acc_1 - b6[1]);

		for (int j = 0; j < 3 * L1_horizon + 2; j++) {
			temp_z_k_1[j] = z_k[j] + (v_k_1[j] - y_k[j]) / rho;
		}

		for (int j = 0; j < L1_horizon; j++) {
			x_k_1[j] 					= alfa * u[j] + (1 - alfa) * x_k[j];
			x_k_1[j + L1_horizon] 		= alfa * t[j] + (1 - alfa) * x_k[j + L1_horizon];

			z_k_1[j] 					= alfa * temp_z_k_1[j] + (1 - alfa) * z_k[j] + y_k[j] / rho;
			if (z_k_1[j] > 0)
				z_k_1[j] = 0;
			else if (z_k_1[j] < -2)
				z_k_1[j] = -2;

			z_k_1[j + L1_horizon] 		= alfa * temp_z_k_1[j + L1_horizon] + (1 - alfa) * z_k[j + L1_horizon] + y_k[j + L1_horizon] / rho;
			if (z_k_1[j + L1_horizon] > 2)
				z_k_1[j + L1_horizon] = 2;
			else if (z_k_1[j + L1_horizon] < 0)
				z_k_1[j + L1_horizon] = 0;

			z_k_1[j + 2 * L1_horizon] 	= alfa * temp_z_k_1[j + 2 * L1_horizon] + (1 - alfa) * z_k[j + 2 * L1_horizon] + y_k[j + 2 * L1_horizon] / rho;
			if (z_k_1[j + 2 * L1_horizon] > 1)
				z_k_1[j + 2 * L1_horizon] = 1;
			else if (z_k_1[j + 2 * L1_horizon] < 0)
				z_k_1[j + 2 * L1_horizon] = 0;
		}

		z_k_1[3 * L1_horizon] = -state_w0;
		z_k_1[3 * L1_horizon + 1] = - L1_horizon * Ka_d * state_w0 - state_q0;

		for (int j = 0; j < 3 * L1_horizon + 2; j++) {
			y_k_1[j] = y_k[j] + rho * (alfa * temp_z_k_1[j] + (1 - alfa) * z_k[j] - z_k_1[j]);
		}

		for (int j = 0; j < 2 * L1_horizon; j++) {
			x_k[j] = x_k_1[j];
		}

		for (int j = 0; j < 3 * L1_horizon + 2; j++) {
			y_k[j] = y_k_1[j];
			z_k[j] = z_k_1[j];
		}
	}

	for (int i = 0; i < L1_horizon; i++) {
		u_opt[i] = x_k[i];
	}
}
