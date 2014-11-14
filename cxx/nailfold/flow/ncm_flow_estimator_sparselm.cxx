#include "ncm_flow_estimator_sparselm.h"

#include <nailfold/flow/ncm_flow_field.h>
#include <nailfold/flow/ncm_flow_cost_function_sparselm.h>

#include <splm.h>

struct info_struct {
  info_struct(
    double* info_vector,
    int n_values);

  double initial_squared_error;
  double final_squared_error;
  double max_Jt_e;
  double dp_magnitude;
  double mu_normalized;
  double n_iterations;
  double termination_reason;
  double n_function_evals;
  double n_jacobian_evals;
  double n_systems_solved;
};

info_struct::info_struct(
    double* info,
    int n_values)
{
  initial_squared_error = info[0];

  //* info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], 
  //* all computed at estimated p.
  final_squared_error = info[1];
  max_Jt_e = info[2];
  dp_magnitude = info[3];
  mu_normalized = info[4];

  n_iterations = info[5];

  //* info[6]=reason for terminating: 1 - stopped by small gradient J^T e
  //*                                 2 - stopped by small dp
  //*                                 3 - stopped by itmax
  //*                                 4 - singular matrix. Restart from current p with increased mu 
  //*                                 5 - too many failed attempts to increase damping. Restart with increased mu
  //*                                 6 - stopped by small ||e||_2
  //*                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. User error
  termination_reason = info[5];
  n_function_evals = info[6];
  n_jacobian_evals = info[7];
  if (n_values > 8)
    n_systems_solved = info[8];
}
 
//
//  Public methods
//
 
//: Estimate the flow field from the input image stack.
void ncm_flow_estimator_sparselm::estimate(
  ncm_flow_field& flow_field,
  bool update)
{
  ncm_flow_cost_function_sparselm func(
    image_stack_,
    warped_image_stack_,
    robustifier_);

  // Input options
  double opts[SPLM_OPTS_SZ];
  opts[0] = SPLM_INIT_MU; 
  opts[1] = 1e-2; 
  opts[2] = 1e-2;
  opts[3] = 1e-2;
  opts[4] = 1e-3; // relevant only if finite difference approximation to Jacobian is used
  opts[5] = SPLM_CHOLMOD;

  // Output information?
  double info[SPLM_INFO_SZ];

  const unsigned n_parameters = flow_field.n_parameters();
  double* parameter_vector = new double[n_parameters];
  flow_field.get_to_vector(parameter_vector);

  int result = sparselm_dercrs(
    &ncm_flow_cost_function_sparselm::f, 
    &ncm_flow_cost_function_sparselm::crs_jacobian,
    parameter_vector, 
    /* observations = */ NULL,
    n_parameters, 
    /* const int nconvars = */ 0, 
    func.n_observations(), 
    func.Jnnz(), func.JtJnnz(),
    /* int itmax = */ 5, 
    opts, info,
    /* void *adata= */ &func);

  info_struct is(info, SPLM_INFO_SZ);

  if (update)
    flow_field.update_from_vector(parameter_vector);
  else
    flow_field.set_from_vector(parameter_vector);

  delete[] parameter_vector;
}

