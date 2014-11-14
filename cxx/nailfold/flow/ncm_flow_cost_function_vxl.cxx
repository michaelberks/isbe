#include "ncm_flow_cost_function_vxl.h"

#include <vil/vil_math.h>
#include <vil/algo/vil_sobel_3x3.h>

#include <nailfold/flow/ncm_flow_null_robustifier.h>
#include <nailfold/flow/ncm_flow_field.h>
 
//
//  Public methods
//
 
//: Constructor.
ncm_flow_cost_function_vxl::ncm_flow_cost_function_vxl(
  unsigned int num_a,
  unsigned int num_params_per_a,
  unsigned int num_b,
  unsigned int num_params_per_b,
  unsigned int num_params_c,
  unsigned int num_residuals_per_e,
  UseGradient g /* = use_gradient */,
  UseWeights w /* = no_weights */,
  vcl_vector< vil_image_view<vxl_byte> > const* image_stack /* = NULL */,
  ncm_flow_robustifier const* robustifier /* = NULL */)
: vnl_sparse_lst_sqr_function(num_a, num_params_per_a, 
                              num_b, num_params_per_b,
                              num_params_c,
                              num_residuals_per_e,
                              g, w),
  robustifier_(robustifier)
{
  differentiate_image_stack(*image_stack);

  // TODO: something with robustifier?
  
}
 
//: Constructor.
ncm_flow_cost_function_vxl::ncm_flow_cost_function_vxl(
  unsigned int num_a,
  unsigned int num_params_per_a,
  unsigned int num_b,
  unsigned int num_params_per_b,
  unsigned int num_params_c,
  const vcl_vector<vcl_vector<bool> >& xmask,
  unsigned int num_residuals_per_e,
  UseGradient g /* = use_gradient */,
  UseWeights w /* = no_weights */,
  vcl_vector< vil_image_view<vxl_byte> > const* image_stack /* = NULL */,
  ncm_flow_robustifier const* robustifier /* = NULL */)
: vnl_sparse_lst_sqr_function(num_a, num_params_per_a, 
                              num_b, num_params_per_b,
                              num_params_c,
                              xmask,
                              num_residuals_per_e,
                              g, w),
  robustifier_(robustifier)
{
  differentiate_image_stack(*image_stack);

  // TODO: something with robustifier?
  
}
 
//: Destructor.
ncm_flow_cost_function_vxl::~ncm_flow_cost_function_vxl()
{
}
 
//: Compute all residuals.
//  Given the parameter vectors a, b, and c, compute the vector of residuals f.
//  f has been sized appropriately before the call.
//  The default implementation computes f by calling fij for each valid
//  pair of i and j.  You do not need to overload this method unless you
//  want to provide a more efficient implementation for your problem.

/* virtual */ 
void ncm_flow_cost_function_vxl::f(
  vnl_vector<double> const& a,
  vnl_vector<double> const& b,
  vnl_vector<double> const& c,
  vnl_vector<double>& f)
{  
  long f_offset;
  
  ncm_flow_field flow_field_to_test(It_[1].ni(), It_[1].nj());

  flow_field_to_test.set_from_vectors(a, b, c);

  vnl_vector<double> f_in;
  f_in.set_size(n_residuals());

  unsigned n_added = 0;
  
  f_offset = 0;
  n_added = brightness_constancy(flow_field_to_test, f_in, f_offset);

  //f_offset = n_added;
  //n_added = spatial_smoothness(flow_field_to_test, f_in, f_offset);

  f = robustifier_->f(f_in);
}
 
/* virtual */ 
void ncm_flow_cost_function_vxl::fij(
  int i, int j,
  vnl_vector<double> const& aij,
  vnl_vector<double> const& bij,
  vnl_vector<double> const& c,
  vnl_vector<double>& fij)
{  
  long f_offset;
  
  ncm_flow_field flow_field_to_test(It_[1].ni(), It_[1].nj());

  flow_field_to_test.set_from_vectors(a, b, c);

  vnl_vector<double> f_in;
  f_in.set_size(n_residuals());

  unsigned n_added = 0;
  
  f_offset = 0;
  n_added = brightness_constancy(flow_field_to_test, f_in, f_offset);

  //f_offset = n_added;
  //n_added = spatial_smoothness(flow_field_to_test, f_in, f_offset);

  f = robustifier_->f(f_in);
}
 
//
//  Private methods
//

void ncm_flow_cost_function_vxl::differentiate_image_stack(
  vcl_vector< vil_image_view<vxl_byte> > const& I)
{
  const unsigned n_images = I.size();

  Ix_.resize(n_images);
  Iy_.resize(n_images);
  It_.resize(n_images);

  const unsigned ni = I[0].ni();
  const unsigned nj = I[0].nj();

  // Differentiate images in x, y and t
  for (unsigned i = 1; i < n_images-1; ++i)
  {
    Ix_[i].set_size(ni, nj);
    Iy_[i].set_size(ni, nj);
    It_[i].set_size(ni, nj);

    // Approximate spatial derivatives.
    vil_sobel_3x3(I[i],
                  Ix_[i], Iy_[i]);

    // Approximate temporal derivatives.
    vil_math_image_difference(I[i+1], I[i-1],
                              It_[i]);
    vil_math_scale_values(It_[i], 0.5);
  }
}
 
//: Compute the brightness constancy penalties
unsigned long ncm_flow_cost_function_vxl::brightness_constancy(
  ncm_flow_field const& flow_field_to_test,
  vnl_vector<double>& f,
  long f_offset)
{
  unsigned n_added = 0;

  const unsigned ni = flow_field_to_test.u().ni();
  const unsigned nj = flow_field_to_test.u().nj();
  const unsigned n_images = It_.size();

  for (unsigned i = 0; i < ni; ++i)
  {
    for (unsigned j = 0; j < nj; ++j)
    {
      for (unsigned k = 1; k < n_images-1; ++k)
      {
        double rx = Ix_[k](i,j) * flow_field_to_test.u(i,j);
        double ry = Iy_[k](i,j) * flow_field_to_test.v(i,j);
        double rt = It_[k](i,j);

        f[f_offset] = (rx + ry) + rt;
        ++f_offset;

        ++n_added;
      }
    }
  }

  return n_added;
}
 
//: Compute the spatial smoothness penalties
unsigned long ncm_flow_cost_function_vxl::spatial_smoothness(
  ncm_flow_field const& flow_field_to_test,
  vnl_vector<double>& f,
  long f_offset)
{
  unsigned n_added = 0;

  const unsigned ni = flow_field_to_test.u().ni();
  const unsigned nj = flow_field_to_test.u().nj();

  for (unsigned i = 0; i < ni; ++i)
  { 
    for (unsigned j = 0; j < nj; ++j)
    {
      if (i > 0)
      {
        f[f_offset] = flow_field_to_test.u(i,j) - flow_field_to_test.u(i-1,j);
        ++f_offset;

        f[f_offset] = flow_field_to_test.v(i,j) - flow_field_to_test.v(i-1,j);
        ++f_offset;

        n_added += 2;
      }

      if (j > 0)
      {
        f[f_offset] = flow_field_to_test.u(i,j) - flow_field_to_test.u(i,j-1);
        ++f_offset;

        f[f_offset] = flow_field_to_test.v(i,j) - flow_field_to_test.v(i,j-1);
        ++f_offset;

        n_added += 2;
      }
    }
  }

  return n_added;
}
 
//: Return number of observations we'll take into account.
unsigned long ncm_flow_cost_function_vxl::n_observations() const
{
  bool use_observation_mask_ = false;

  if (use_observation_mask_)
  {
    unsigned long observation_mask_count_ = 0;
    return observation_mask_count_;
  }
  else
  {
    const unsigned n_images = It_.size();
  
    if (n_images > 2)
    {
      const unsigned ni = It_[1].ni();
      const unsigned nj = It_[1].nj();

      return (ni) * (nj) * (n_images-2);
    }
    else
      return 0;
  }
}
 
//: Return the total number of residuals to optimize.
unsigned long ncm_flow_cost_function_vxl::n_residuals() const
{
  const unsigned ni = It_[1].ni();
  const unsigned nj = It_[1].nj();
  const unsigned n_images = It_.size();

  const unsigned long n_brightness = n_observations();

  //const unsigned long n_smoothness = (ni-1) *  nj    * 2 + 
  //                                    ni    * (nj-1) * 2;
  const unsigned long n_smoothness = 0;

  return n_brightness + n_smoothness;
}