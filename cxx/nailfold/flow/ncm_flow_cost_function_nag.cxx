#include "ncm_flow_cost_function_nag.h"

#include <vil/vil_math.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_convolve_1d.h>

#include <nailfold/flow/ncm_flow_null_robustifier.h>
#include <nailfold/flow/ncm_flow_field.h>
 
//
//  Public methods
//
 
//: Constructor.
ncm_flow_cost_function_nag::ncm_flow_cost_function_nag(
  vcl_vector< vil_image_view<float> > const* image_stack,
  vcl_vector< vil_image_view<float> > const* warped_image_stack /* = NULL */,
  ncm_flow_robustifier const* robustifier /* = NULL */)
: ncm_flow_cost_function(image_stack, warped_image_stack, robustifier)
{
  // All of the construction is handled by the base class.
}
 
//: Destructor.
ncm_flow_cost_function_nag::~ncm_flow_cost_function_nag()
{
}
 
