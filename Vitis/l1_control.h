#define Ts 6.283185307179586

//#define max_size 300
const int max_size = 300;

#define Ka 1
#define Kb 0.0014

#define Ka_d 6.283185307179586
#define Kb1_d 0.008796459430051
#define Kb2_d 0.027634892323050

#define Kb1_d_sq (Kb1_d * Kb1_d)
#define Kb2_d_sq (Kb2_d * Kb2_d)

#define inv_Ka_Kb (1.0 / (Ka * Kb))
#define inv_Kb (1.0 / Kb)

#define rho 3
#define sigma 0.01
#define alfa 1.95

#define inv_sigma_3_rho (1.0 / (sigma + 3 * rho))

#define constant_diag_value (sqrt(2.0 + sigma / rho))

#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <ap_axi_sdata.h>
#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_math.h>
//#include "polinearsolver.h"
//#include "polinearsolver_orig.h"
#include "chol_decomp.h"
#include "ltris.h"

typedef ap_axis<64, 0, 0, 0> io_variable_axis_wrapper;

#define io_variable_size 64
#define io_variable_integer_size 32

typedef ap_fixed<io_variable_size, io_variable_integer_size> io_variable;

#define state_variable_size 55
#define state_variable_integer_size 27

typedef ap_fixed<state_variable_size, state_variable_integer_size> state_variable;

#define ADMM_variable_size 40
#define ADMM_variable_integer_size 20

typedef ap_fixed<ADMM_variable_size, ADMM_variable_integer_size> ADMM_variable;

#define constant_variable_size 32
#define constant_variable_integer_size 16

typedef ap_fixed<constant_variable_size, constant_variable_integer_size> constant_variable;

#define inv_rho (1.0 / rho)

#define r 0.5

#define n_of_test_cases 1

#define n_of_iterations 700

void hoc_module(hls::stream<io_variable_axis_wrapper> &in_state, hls::stream<io_variable_axis_wrapper> &out_port);
//void hoc_module(io_variable *in_state, io_variable *out_port);
