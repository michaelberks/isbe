#include "ncm_flow_estimator_nag.h"

#include <nailfold/flow/ncm_flow_field.h>
#include <nailfold/flow/ncm_flow_cost_function_nag.h>
 
//
//  Public methods
//
 
//: Estimate the flow field from the input image stack.
void ncm_flow_estimator_nag::estimate(
  ncm_flow_field& flow_field,
  bool update)
{
  ncm_flow_cost_function_nag func(
    image_stack_,
    warped_image_stack_,
    robustifier_);

  const unsigned n_parameters = flow_field.n_parameters();
  double* parameter_vector = new double[n_parameters];
  flow_field.get_to_vector(parameter_vector);

  // Determine the optimal parameters

  if (update)
    flow_field.update_from_vector(parameter_vector);
  else
    flow_field.set_from_vector(parameter_vector);

  delete[] parameter_vector;
}

