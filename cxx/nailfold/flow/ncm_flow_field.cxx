#include "ncm_flow_field.h"

#include <vcl_cmath.h>

#include <vil/vil_math.h>
#include <vil/vil_convert.h>

#include <vil/algo/vil_colour_space.h>
 
//:
ncm_flow_field::ncm_flow_field(
  unsigned ni /* = 0 */,
  unsigned nj /* = 0 */)
{
  set_size(ni, nj);
}

//:
ncm_flow_field::~ncm_flow_field()
{
}
 
//:
void ncm_flow_field::set_size(
  unsigned ni, 
  unsigned nj)
{
  u_.set_size(ni, nj);
  u_.fill(0.0);

  v_.set_size(ni, nj);
  v_.fill(0.0);

  //covariance_.set_size(ni, nj);
  //covariance_.fill(0.0);
}

//: Set contraints/priors on the orientation field.
void ncm_flow_field::set_orientation()
{
}
 
//: Set constraints/priors over the vesselness (or presence) field.
void ncm_flow_field::set_presence()
{
}
 
//: Set constraints/priors over the interframe displacements (i.e. jitter).
void ncm_flow_field::set_displacements()
{
}
  
//: Return the horizontal flow field.
vil_image_view<double> const& ncm_flow_field::u() const
{
  return u_;
}
 
//: Return the vertical flow field.
vil_image_view<double> const& ncm_flow_field::v() const
{
  return v_;
}
 
//: Return the horizontal flow at a given location.
double ncm_flow_field::u(unsigned i, unsigned j) const
{
  return u_(i,j);
}
 
//: Return the vertical flow at a given location.
double ncm_flow_field::v(unsigned i, unsigned j) const
{
  return v_(i,j);
}
 
//: Return the index in the parameter vector of the horizontal flow at (i,j).
int ncm_flow_field::index_of_u(unsigned i, unsigned j)
{
  return 2 * (j * u_.ni() + i);
}
 
//: Return the index in the parameter vector of the vertical flow at (i,j).
int ncm_flow_field::index_of_v(unsigned i, unsigned j)
{
  return 2 * (j * u_.ni() + i) + 1;
}
 
//:
vil_image_view<double> ncm_flow_field::magnitude()
{
  return vil_image_view<double>(0,0);
}
 
//:
vil_image_view<double> ncm_flow_field::direction()
{
  return vil_image_view<double>(0,0);
}
 
//:
vil_image_view<vxl_byte> ncm_flow_field::as_colormap()
{
  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();

  double max_magnitude_cached = max_magnitude();

  const double bg_val = 1.0;

  vil_image_view<double> rgbd(ni, nj, 3);

  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      double angle = vcl_atan2(-v_(i,j), -u_(i,j));
      double magnitude = vcl_sqrt(u_(i,j)*u_(i,j) + v_(i,j)*v_(i,j)) / 
                         max_magnitude_cached;

      double r, g, b;
      vil_colour_space_HSV_to_RGB(
        angle*(180.0/3.141592), 
        (1.0-bg_val) + bg_val*magnitude, 
        255.0 * (bg_val + (1.0 - bg_val)*magnitude),

        &r, &g, &b);

      rgbd(i,j,0) = r;
      rgbd(i,j,1) = g;
      rgbd(i,j,2) = b;
    }
  }

  vil_image_view<vxl_byte> rgb(ni, nj, 3);
  vil_convert_cast(rgbd, rgb);

  return rgb;
}
 
//: Return the total number of parameters.
unsigned ncm_flow_field::n_parameters()
{
  return u_.ni() * u_.nj() + 
         v_.ni() * v_.nj();
}

//: Copy the flow field parameters from three vectors a, b and c (used in
//  nonlinear optimization).
void ncm_flow_field::set_from_vector(
  double const* v)
{
  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();

  unsigned v_index = 0;
  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      u_(i,j) = v[v_index];
      ++v_index;

      v_(i,j) = v[v_index];
      ++v_index;
    }
  }
}
 
//: Copy the flow field parameters from three vectors a, b and c (used in
//  nonlinear optimization).
void ncm_flow_field::update_from_vector(
  double* v)
{
  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();

  unsigned v_index = 0;
  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      u_(i,j) += v[v_index];
      ++v_index;

      v_(i,j) += v[v_index];
      ++v_index;
    }
  }
}
 
//: Copy the flow field parameters to three vectors a, b and c (used in
//  nonlinear optimization).
void ncm_flow_field::get_to_vector(
  double* v)
{
  //v.set_size(n_parameters());

  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();
  
  unsigned v_index = 0;
  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      v[v_index] = u_(i,j);
      ++v_index;

      v[v_index] = v_(i,j);
      ++v_index;
    }
  }
}

//: Upscale the flow map by a given interpolation factor.
//  Values are replicated rather than interpolated. This won't matter much in
//  practice as the upscaling is used to initialize another round of nonlinear
//  optimization.
void ncm_flow_field::upscale(unsigned factor /* = 2 */)
{
  // Do nothing if too small a factor is specified.
  if (factor < 2)
    return;

  // Upscale the flow field (spatial flow)
  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();

  vil_image_view<double> new_u(ni*factor, nj*factor);
  vil_image_view<double> new_v(ni*factor, nj*factor);
  
  for (unsigned j = 0; j < nj*factor; ++j)
  {
    for (unsigned i = 0; i < ni*factor; ++i)
    {
      new_u(i,j) = factor * u_(i/factor, j/factor);
      new_v(i,j) = factor * v_(i/factor, j/factor);
    }
  }

  u_ = new_u;
  v_ = new_v;

  // Upscale the temporal flow (including the displacements)
}
 
//
//  Private methods
//
 
//: Return the greatest flow magnitude across the field
double ncm_flow_field::max_magnitude()
{
  const unsigned ni = u_.ni();
  const unsigned nj = u_.nj();

  double max_magnitude = 0.0;

  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      const double magnitude_ij = 
        vcl_sqrt(u_(i,j)*u_(i,j) + v_(i,j)*v_(i,j));

      if (magnitude_ij > max_magnitude)
        max_magnitude = magnitude_ij;
    }
  }

  return max_magnitude;
}
