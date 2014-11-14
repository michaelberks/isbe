#include "ncm_flow_cost_function.h"

#include <vil/vil_math.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_convolve_1d.h>

#include <nailfold/flow/ncm_flow_null_robustifier.h>
#include <nailfold/flow/ncm_flow_field.h>
 
//
//  Public methods
//
 
//: Constructor.
ncm_flow_cost_function::ncm_flow_cost_function(
  vcl_vector< vil_image_view<float> > const* image_stack,
  vcl_vector< vil_image_view<float> > const* warped_image_stack /* = NULL */,
  ncm_flow_robustifier const* robustifier /* = NULL */)
: alpha_(0.8),
  robustifier_(robustifier)
{
  differentiate_spatially(*image_stack);

  if (warped_image_stack == NULL)
    differentiate_temporally(*image_stack, *image_stack);
  else
    differentiate_temporally(*image_stack, *warped_image_stack);

  // Cache the image stack size.
  ni_ = It_[1].ni();
  nj_ = It_[1].nj();
  n_images_ = It_.size();

  // TODO: something with robustifier?
  
}
 
//: Destructor.
ncm_flow_cost_function::~ncm_flow_cost_function()
{
}
 
//: Return the number of parameters (variables).
int ncm_flow_cost_function::n_parameters()
{
  return -1;
}
 
//: Return all parameters as a vector of doubles.
double* ncm_flow_cost_function::parameters()
{
  //vnl_vector<double> parameters;

  //flow_field.get_to_vector(parameters);

  //return parameters.data_block();

  return NULL;
}
 
//: Return the number of parameters whose value is fixed.
int ncm_flow_cost_function::n_constrained_parameters()
{
  return 0;
}
 
//: Return the number of observations.
int ncm_flow_cost_function::n_observations()
{
  return n_obs_brightness() + 
         n_obs_smoothness();
}
 
//: Return all observations as a vector of doubles.
double* ncm_flow_cost_function::observations()
{
  return NULL;
}
 
//: Return the number of nonzero elements in the Jacobian.
int ncm_flow_cost_function::Jnnz()
{
  int Jnnz = 0;

  Jnnz += n_obs_brightness() * 2;
  Jnnz += n_obs_smoothness() * 2;

  return Jnnz;
}
 
//: Return the number of nonzero elements in the approximation to the Hessian.
int ncm_flow_cost_function::JtJnnz()
{
  return -1;
}
 
//
//  Private methods
//
 
//:
void ncm_flow_cost_function::smoothe_x(
  vil_image_view<float> const& in,
  vil_image_view<float>& out)
{
  float kernel[3] = {0.25f, 0.5f, 0.25f};

  float const* in_row_ptr = in.top_left_ptr();
  float* out_row_ptr = out.top_left_ptr();

  for (unsigned j = 0; j < in.nj(); ++j)
  {
    float acc = 0.0f;

    vil_convolve_1d(
      in_row_ptr, in.ni(), in.istep(),
      out_row_ptr, out.istep(),
      &kernel[1],
      -1, 1,
      acc,
      vil_convolve_constant_extend,
      vil_convolve_constant_extend);

    in_row_ptr += in.jstep();
    out_row_ptr += out.jstep();
  }
}
 
//:
void ncm_flow_cost_function::smoothe_y(
  vil_image_view<float> const& in,
  vil_image_view<float>& out)
{
  float kernel[3] = {0.25f, 0.5f, 0.25f};

  float const* in_col_ptr = in.top_left_ptr();
  float* out_col_ptr = out.top_left_ptr();

  for (unsigned i = 0; i < in.ni(); ++i)
  {
    float acc = 0.0f;

    vil_convolve_1d(
      in_col_ptr, in.nj(), in.jstep(),
      out_col_ptr, out.jstep(),
      &kernel[1],
      -1, 1,
      acc,
      vil_convolve_constant_extend,
      vil_convolve_constant_extend);

    in_col_ptr += in.istep();
    out_col_ptr += out.istep();
  }
}
 
//:
void ncm_flow_cost_function::differentiate_spatially(
  vcl_vector< vil_image_view<float> > const& I)
{
  const unsigned n_images = I.size();

  Ix_.resize(n_images);
  Iy_.resize(n_images);

  const unsigned ni = I[0].ni();
  const unsigned nj = I[0].nj();

  vil_image_view<float> tmp_img(ni, nj);

  for (unsigned i = 0; i < n_images; ++i)
  {
    Ix_[i].set_size(ni, nj);
    Iy_[i].set_size(ni, nj);

    // Approximate spatial derivatives.
    vil_sobel_3x3(I[i],
                  Ix_[i], Iy_[i]);

    smoothe_x(Ix_[i], tmp_img);
    Ix_[i].deep_copy(tmp_img);

    smoothe_y(Iy_[i], tmp_img);
    Iy_[i].deep_copy(tmp_img);
  }
}
 
//: Compute temporal derivatives between I0(f+1) and I(f).
//  Passing the same sequence for both I0 and I1 gives the usual derivatives
//  as computed from consecutive frames.
//  Warping a sequence (e.g. according to an estimated flow field) permits
//  multiresolution flow estimation.
void ncm_flow_cost_function::differentiate_temporally(
  vcl_vector< vil_image_view<float> > const& I0,
  vcl_vector< vil_image_view<float> > const& I1)
{
  const unsigned n_images = I0.size();

  It_.resize(n_images);

  const unsigned ni = I0[0].ni();
  const unsigned nj = I0[0].nj();

  vil_image_view<float> tmp_img(ni, nj);

  for (unsigned i = 0; i < n_images-1; ++i)
  {
    vil_math_image_difference(I1[i+1], I0[i],
                              It_[i]);

    // Smoothe spatially - questionable whether this is necessary
    // (Spatial derivatives aren't smoothed temporally.)
    smoothe_x(It_[i], tmp_img);
    smoothe_y(tmp_img, It_[i]);
  }

  // Replicate finite differencing for final frame.
  It_[n_images-1].deep_copy(It_[n_images-2]);

  // Smoothe frames 1..n_images-2 temporally to give something 
  // roughly equivalent to (I[f+1]-I[f-1])/2 in practice.
  // Note, this is done in reverse to avoid recursive averaging.
  for (unsigned i = n_images-2; i > 0; --i)
  {
    vil_math_image_sum(It_[i], It_[i-1], It_[i]);
    vil_math_scale_values(It_[i], 0.5);
  }
}
 
//: Compute the brightness constancy penalties
unsigned long ncm_flow_cost_function::brightness_constancy(
  ncm_flow_field const& flow_field_to_test,
  double* f,
  long f_offset0)
{
  long f_offset = f_offset0;
  unsigned n_added = 0;

  const double weight = 0.5 * alpha_;

  for (unsigned k = 0; k < n_images_; ++k)
  {
    vil_image_view<float> Ix_k = Ix_[k];
    vil_image_view<float> Iy_k = Iy_[k];
    vil_image_view<float> It_k = It_[k];

    for (unsigned j = 0; j < nj_; ++j)
    {
      for (unsigned i = 0; i < ni_; ++i)
      {
        const double rx = Ix_k(i,j) * flow_field_to_test.u(i,j);
        const double ry = Iy_k(i,j) * flow_field_to_test.v(i,j);
        const double rt = It_k(i,j);

        f[f_offset] = (rx + ry) + rt;
        f[f_offset] *= weight;

        f_offset += n_images_;

        ++n_added;
      }
    }

    ++f_offset0;
    f_offset = f_offset0;
  }

  return n_added;
}
 
//: Compute the spatial smoothness penalties
unsigned long ncm_flow_cost_function::spatial_smoothness(
  ncm_flow_field const& flow_field_to_test,
  double* f,
  long f_offset)
{
  unsigned n_added = 0;

  const double weight_hv = (1.0 - alpha_) * 
                           static_cast<double>(n_images_) / 6.0;

  for (unsigned j = 0; j < nj_; ++j)
  {
    for (unsigned i = 0; i < ni_; ++i)
    { 
      if (i < ni_-1)
      {
        f[f_offset] = flow_field_to_test.u(i,j) - flow_field_to_test.u(i+1,j);
        f[f_offset] *= weight_hv;
        ++f_offset;

        f[f_offset] = flow_field_to_test.v(i,j) - flow_field_to_test.v(i+1,j);
        f[f_offset] *= weight_hv;
        ++f_offset;

        n_added += 2;
      }

      if (j < nj_-1)
      {
        f[f_offset] = flow_field_to_test.u(i,j) - flow_field_to_test.u(i,j+1);
        f[f_offset] *= weight_hv;
        ++f_offset;

        f[f_offset] = flow_field_to_test.v(i,j) - flow_field_to_test.v(i,j+1);
        f[f_offset] *= weight_hv;
        ++f_offset;

        n_added += 2;
      }
    }
  }

  return n_added;
}
 
//: Return number of brightness constraints.
unsigned long ncm_flow_cost_function::n_obs_brightness() const
{
  bool use_observation_mask_ = false;

  if (use_observation_mask_)
  {
    unsigned long observation_mask_count_ = 0;
    return observation_mask_count_;
  }
  else
  { 
    return ni_ * nj_ * n_images_;
  }
}
 
//: Return the number of smoothness constraints.
unsigned long ncm_flow_cost_function::n_obs_smoothness() const
{
  const unsigned n_horizontal = (ni_-1) * nj_;
  const unsigned n_vertical = ni_ * (nj_-1);

  // One horizontal and vertical for u and the same for v.
  return 2 * (n_horizontal + n_vertical);
}
