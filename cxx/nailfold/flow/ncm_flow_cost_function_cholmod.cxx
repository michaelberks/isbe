#include "ncm_flow_cost_function_cholmod.h"

#include <vil/vil_math.h>

#include <nailfold/flow/ncm_flow_null_robustifier.h>
#include <nailfold/flow/ncm_flow_field.h>
 
//
//  Public methods
//
 
//: Constructor.
ncm_flow_cost_function_cholmod::ncm_flow_cost_function_cholmod(
    vcl_vector< vil_image_view<float> > const* image_stack,
    vcl_vector< vil_image_view<float> > const* warped_image_stack /* = NULL */,
    ncm_flow_robustifier const* robustifier /* = NULL */)
: ncm_flow_cost_function(image_stack, warped_image_stack, robustifier)
{
}
 
//: Destructor.
ncm_flow_cost_function_cholmod::~ncm_flow_cost_function_cholmod()
{
}
 
vnl_vector<double> ncm_flow_cost_function_cholmod::f(
  vnl_vector<double> const& parameters)
{
  ncm_flow_field flow_field_to_test(ni_, nj_);

  flow_field_to_test.set_from_vector(parameters.data_block());

  vnl_vector<double> x_hat(n_observations());

  unsigned n_added = 0;
  long f_offset = 0;

  n_added = brightness_constancy(flow_field_to_test, x_hat.data_block(), f_offset);
  f_offset += n_added;
  
  n_added = spatial_smoothness(flow_field_to_test, x_hat.data_block(), f_offset);
  f_offset += n_added;

  return x_hat;
}
 
vnl_sparse_matrix<double> ncm_flow_cost_function_cholmod::jacobian(
  vnl_vector<double> const& parameters)
{
  vnl_sparse_matrix<double> J(n_observations(), 
                              parameters.size());

  // Add the brightness gradients.
  unsigned u_index = 0;
  unsigned v_index = 1;
  unsigned row = 0;

  double weight = 0.5 * alpha_;

  for (unsigned j = 0; j < nj_; ++j)
  {
    for (unsigned i = 0; i < ni_; ++i)
    {
      for (unsigned k = 0; k < n_images_; ++k)
      {
        J(row, u_index) = weight * Ix_[k](i,j);
        J(row, v_index) = weight * Iy_[k](i,j);

        ++row;
      }

      u_index += 2;
      v_index += 2;
    }
  }

  // Add the smoothness gradients.
  ncm_flow_field ff(ni_, nj_);
  const double weight_hv = (1.0 - alpha_) * 
                           static_cast<double>(n_images_) / 6.0;

  for (unsigned j = 0; j < nj_; ++j)
  {
    for (unsigned i = 0; i < ni_; ++i)
    {
      if (i < ni_-1)
      {
        J(row, ff.index_of_u(i,j  )) =  weight_hv;
        J(row, ff.index_of_u(i+1,j)) = -weight_hv;
        ++row;

        J(row, ff.index_of_v(i,j  )) =  weight_hv;
        J(row, ff.index_of_v(i+1,j)) = -weight_hv;
        ++row;
      }

      if (j < nj_-1)
      {
        J(row, ff.index_of_u(i,j  )) =  weight_hv;
        J(row, ff.index_of_u(i,j+1)) = -weight_hv;
        ++row;

        J(row, ff.index_of_v(i,j  )) =  weight_hv;
        J(row, ff.index_of_v(i,j+1)) = -weight_hv;
        ++row;
      }
    }
  }

  return J;
}