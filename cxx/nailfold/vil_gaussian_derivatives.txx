#ifndef __vil_gaussian_derivatives_txx__
#define __vil_gaussian_derivatives_txx__
//:
// \file
// \brief Filter images with first (second) derivative of a Gaussian, and
//        return the edge (line) strength and orientation at every pixel

// \author Phil Tresadern

#include "vil_gaussian_derivatives.h"

#include <vcl_cmath.h>

#include <vil/vil_image_view.h>
#include <vil/vil_transform.h>
#include <vil/vil_math.h>
#include <vil/algo/vil_gauss_filter.h>

// functors for common trig functions

class vil_math_atan_functor
{
 public:
  vxl_byte operator()(vxl_byte x) const { return static_cast<vxl_byte>(0.5+vcl_atan(double(x))); }
  unsigned operator()(unsigned x) const { return static_cast<unsigned int>(0.5+vcl_atan(double(x))); }
  int operator()(int x)           const { return static_cast<int>(0.5+vcl_atan(double(x))); }
  short operator()(short x)       const { return static_cast<short>(0.5+vcl_atan(double(x))); }
  float operator()(float x)       const { return vcl_atan(x); }
  double operator()(double x)     const { return vcl_atan(x); }
};
class vil_math_atan2_functor
{
 public:
  vxl_byte operator()(vxl_byte y, vxl_byte x) const 
    { return static_cast<vxl_byte>(0.5+vcl_atan2(double(y),double(x))); }
  unsigned operator()(unsigned y, unsigned x) const 
    { return static_cast<unsigned int>(0.5+vcl_atan2(double(y),double(x))); }
  int operator()(int y, int x) const 
    { return static_cast<int>(0.5+vcl_atan2(double(y),double(x))); }
  short operator()(short y, short x) const 
    { return static_cast<short>(0.5+vcl_atan2(double(y),double(x))); }
  float operator()(float y, float x) const 
    { return vcl_atan2(y,x); }
  double operator()(double y, double x) const 
    { return vcl_atan2(y,x); }
};
class vil_math_cos_functor
{
 public:
  vxl_byte operator()(vxl_byte x) const { return static_cast<vxl_byte>(0.5+vcl_cos(double(x))); }
  unsigned operator()(unsigned x) const { return static_cast<unsigned int>(0.5+vcl_cos(double(x))); }
  int operator()(int x)           const { return static_cast<int>(0.5+vcl_cos(double(x))); }
  short operator()(short x)       const { return static_cast<short>(0.5+vcl_cos(double(x))); }
  float operator()(float x)       const { return vcl_cos(x); }
  double operator()(double x)     const { return vcl_cos(x); }
};
class vil_math_sin_functor
{
 public:
  vxl_byte operator()(vxl_byte x) const { return static_cast<vxl_byte>(0.5+vcl_sin(double(x))); }
  unsigned operator()(unsigned x) const { return static_cast<unsigned int>(0.5+vcl_sin(double(x))); }
  int operator()(int x)           const { return static_cast<int>(0.5+vcl_sin(double(x))); }
  short operator()(short x)       const { return static_cast<short>(0.5+vcl_sin(double(x))); }
  float operator()(float x)       const { return vcl_sin(x); }
  double operator()(double x)     const { return vcl_sin(x); }
};

// First derivatives

//
//: Compute gradient strength and orientation from first derivatives of a
//  Gaussian
template <class srcT, class destT>
void vil_gaussian_1st_derivative(const vil_image_view<srcT>& src,
                                 vil_image_view<destT>& line_strength,
                                 vil_image_view<destT>& orientation,
                                 double sd, /* = 3.0 */ 
                                 int half_width /* = -1*/,
                                 vil_convolve_boundary_option bo /* = vil_convolve_zero_extend */)
{
  const destT half_pi = static_cast<destT>(1.5707963267);
  const destT two_pi = static_cast<destT>(6.2831853071);

  // if half_width is negative then set it to 5*sd
  if (half_width < 0)
    half_width = static_cast<int>(5.0 * sd);

  // compute gaussian kernels
  vcl_vector<double> g(2*half_width+1);
  vil_gauss_filter_gen_ntap(sd, 0, g);

  vcl_vector<double> dg(2*half_width+1);
  vil_gauss_filter_gen_ntap(sd, 1, dg);
  
  // set outputs to correct size
  line_strength.set_size(src.ni(), src.nj());
  orientation.set_size(src.ni(), src.nj());

  // create workspace image for reuse throughout
  vil_image_view<destT> work_im(src.ni(), src.nj(), src.nplanes());
  vil_image_view<destT> work_im_t = vil_transpose(work_im);

  // response to gradient operator in x
  vil_image_view<destT> Ix(src.ni(), src.nj(), src.nplanes());
  {
    vil_convolve_1d(src,work_im,
                    &dg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
    vil_image_view<destT> dest_t = vil_transpose(Ix);
    vil_convolve_1d(work_im_t,dest_t,
                    &g[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
  }

  // response to gradient operator in y
  vil_image_view<destT> Iy(src.ni(), src.nj(), src.nplanes());
  {
    vil_convolve_1d(src,work_im,
                    &g[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
    vil_image_view<destT> dest_t = vil_transpose(Iy);
    vil_convolve_1d(work_im_t,dest_t,
                    &dg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
  }
  vil_math_scale_values(Iy, -1.0);

  // compute angle of greatest response
  //   theta = atan2(Iy, Ix);
  vil_transform(Iy, Ix, orientation, vil_math_atan2_functor());

  // compute strength of response
  // I use references here to make it clearer what is going on, even though the
  // operations are done in-place
  vil_image_view<destT>& Ix_squared = Ix;
  vil_image_view<destT>& Iy_squared = Iy;
  vil_math_image_product(Ix, Ix, Ix_squared);
  vil_math_image_product(Iy, Iy, Iy_squared);
  
  vil_image_view<destT>& line_strength_squared = line_strength;
  vil_math_image_sum(Ix_squared, Iy_squared, line_strength_squared);
  vil_math_sqrt(line_strength_squared);
}                          

// Second derivatives

template <class destT>
static 
void w_response(const vil_image_view<destT>& Ixx,
                const vil_image_view<destT>& Iyy,
                const vil_image_view<destT>& Ixy,
                const vil_image_view<destT>& theta,
                vil_image_view<destT>& resp_parallel,
                vil_image_view<destT>& resp_perpendicular)
{
  //cc = cos(theta_a).^2;
  vil_image_view<destT> cc; 
  cc.deep_copy(theta);
  vil_transform(cc,vil_math_cos_functor());
  vil_math_image_product(cc, cc, cc);

  //ss = sin(theta_a).^2;
  // applying the rule sin*sin = 1-cos*cos is quicker than recomputing
  // sin*sin from trigonometry
  vil_image_view<destT> ss; 
  ss.deep_copy(cc);
  vil_math_scale_and_offset_values(ss,-1,1);

  //s2 = sin(2*theta_a);
  vil_image_view<destT> s2; 
  s2.deep_copy(theta);
  vil_math_scale_values(s2, 2.0);
  vil_transform(s2,vil_math_sin_functor());

  vil_image_view<destT> work_im;

  //wo_theta_a = Ixx.*cc + Iyy.*ss + Ixy.*s2;
  resp_parallel.set_size(theta.ni(), theta.nj());
  vil_math_image_product(Ixx,cc,resp_parallel);
  vil_math_image_product(Iyy,ss,work_im);
  vil_math_image_sum(resp_parallel,work_im,resp_parallel);
  vil_math_image_product(Ixy,s2,work_im);
  vil_math_image_sum(resp_parallel,work_im,resp_parallel);

  // cos(t + pi/2) = -sin(t) => cc(t + pi/2) = ss
  // sin(t + pi/2) = cos(t)  => ss(t + pi/2) = cc
  // sin(2(t + pi/2)) = sin(2t + pi) = -sin(2t) = -s2

  //wo_theta_b = Ixx.*ss + Iyy.*cc - Ixy.*s2;
  resp_perpendicular.set_size(theta.ni(), theta.nj());
  vil_math_image_product(Ixx,ss,resp_perpendicular);
  vil_math_image_product(Iyy,cc,work_im);
  vil_math_image_sum(resp_perpendicular,work_im,resp_perpendicular);
  vil_math_image_product(Ixy,s2,work_im);
  vil_math_image_difference(resp_perpendicular,work_im,resp_perpendicular);
}

//
//:
template <class srcT, class destT>
void vil_gaussian_2nd_derivative(const vil_image_view<srcT>& src,
                                 vil_image_view<destT>& line_strength,
                                 vil_image_view<destT>& orientation,
                                 double sd, /* = 3.0 */ 
                                 int half_width /* = -1*/,
                                 vil_convolve_boundary_option bo /* = vil_convolve_zero_extend */)
{
  const destT half_pi = static_cast<destT>(1.5707963267);
  const destT two_pi = static_cast<destT>(6.2831853071);

  // if half_width is negative then set it to 5*sd
  if (half_width < 0)
    half_width = static_cast<int>(5.0 * sd);
    
  const int full_width = 2*half_width + 1;

  vcl_vector<double> g(full_width);
  vcl_vector<double> dg(full_width);
  vcl_vector<double> ddg(full_width);
  
/*
  // These functions scale the kernels in a funny way, making them incorrect
  // for the oriented filtering.
  
  // compute gaussian kernels
  vil_gauss_filter_gen_ntap(sd, 0, g);
  vil_gauss_filter_gen_ntap(sd, 1, dg);
  vil_gauss_filter_gen_ntap(sd, 2, ddg);
*/

  // Compute gaussian kernels from first principles.
  double k = vcl_pow(static_cast<double>(two_pi), -0.5);
  double sigmasq = sd*sd;
  for (int x = -half_width; x <= half_width; ++x)
  {
    unsigned v = x + half_width;
    g[v]   = k * vcl_exp(-0.5 * (x*x) / sigmasq);
    dg[v]  =  -x / sigmasq * g[v];
    ddg[v] = (-1 / sigmasq * g[v]) - (x / sigmasq * dg[v]);
  }
  
  // set outputs to correct size
  line_strength.set_size(src.ni(), src.nj());
  orientation.set_size(src.ni(), src.nj());

  // create workspace image for reuse throughout
  vil_image_view<destT> work_im(src.ni(), src.nj(), src.nplanes());
  vil_image_view<destT> work_im_t = vil_transpose(work_im);

  // filter source image three times
  vil_image_view<destT> Ixx(src.ni(), src.nj(), src.nplanes());
  {
    vil_convolve_1d(src,work_im,
                    &ddg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
    vil_image_view<destT> dest_t = vil_transpose(Ixx);
    vil_convolve_1d(work_im_t,dest_t,
                    &g[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
  }

  vil_image_view<destT> Iyy(src.ni(), src.nj(), src.nplanes());
  {
    vil_convolve_1d(src,work_im,
                    &g[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
    vil_image_view<destT> dest_t = vil_transpose(Iyy);
    vil_convolve_1d(work_im_t,dest_t,
                    &ddg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
  }

  vil_image_view<destT> Ixy(src.ni(), src.nj(), src.nplanes());
  {
    vil_convolve_1d(src,work_im,
                    &dg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
    vil_image_view<destT> dest_t = vil_transpose(Ixy);
    vil_convolve_1d(work_im_t,dest_t,
                    &dg[half_width],-int(half_width),half_width,
                    destT(), bo, bo);
  }
  // Invert response later
  //vil_math_scale_values(Ixy, -1.0);


  // compute angle of greatest response
  //   theta = atan(2*Ixy ./ (Ixx-Iyy)) / 2;
  
  vil_math_scale_values(Ixy, -2.0); // Response inverted here
  vil_image_view<destT> Idiff;
  vil_math_image_difference(Ixx, Iyy, Idiff);
  vil_image_view<destT> theta;
  vil_transform(Ixy, Idiff, theta, vil_math_atan2_functor());
  vil_math_scale_values(theta, 0.5);

  // compute maximum responses
  // use trigonometric identities for efficiency
  vil_image_view<destT> resp_parallel(src.ni(), src.nj(), src.nplanes());
  vil_image_view<destT> resp_perpendicular(src.ni(), src.nj(), src.nplanes());
  w_response(Ixx, Iyy, Ixy, theta, 
             resp_parallel, resp_perpendicular);
  
  // if wo_theta_a at this scale so far is biggest, swap it into the main
  // line_strength, orientation and scale matrices - note the theta are
  // normal to the line direction, so for theta_a we swap in theta_b and
  // vice-versa

  // line_strength must be initialized (due to the comparison with 
  // *strength_ptr) whereas orientation need not be
  // TODO: could actually assign with one of the responses if we are considering
  // only a single scale at once
  line_strength.fill(0.0f);

  // the following is bigger than it needs to be on account of being done
  // with pointers
  const destT* para_row_ptr = resp_parallel.top_left_ptr();
  const destT* perp_row_ptr = resp_perpendicular.top_left_ptr();
  const destT* theta_row_ptr = theta.top_left_ptr();
  destT* strength_row_ptr = line_strength.top_left_ptr();
  destT* orientation_row_ptr = orientation.top_left_ptr();

  // cache these values to save repeated function calls
  const unsigned src_ni = src.ni();
  const unsigned src_nj = src.nj();
  const int para_istep = resp_parallel.istep();
  const int perp_istep = resp_perpendicular.istep();
  const int theta_istep = theta.istep();
  const int strength_istep = line_strength.istep();
  const int orientation_istep = orientation.istep();
  const int para_jstep = resp_parallel.jstep();
  const int perp_jstep = resp_perpendicular.jstep();
  const int theta_jstep = theta.jstep();
  const int strength_jstep = line_strength.jstep();
  const int orientation_jstep = orientation.jstep();
  
  for (unsigned j = 0; j < src_nj; ++j)
  {
    const destT* para_ptr = para_row_ptr;
    const destT* perp_ptr = perp_row_ptr;
    const destT* theta_ptr = theta_row_ptr;
    destT* strength_ptr = strength_row_ptr;
    destT* orientation_ptr = orientation_row_ptr;

    for (unsigned i = 0; i < src_ni; ++i)
    {
	  destT abs_strength = vcl_abs(*strength_ptr);
	  
	  // TODO: rewrite this to avoid updating variables twice unnecessarily
	  // i.e. if (b > a) { if (c > b) { assign to c, else assign to b }, else assign to a }
      if (vcl_abs(*para_ptr) > abs_strength)
      {
        *strength_ptr = *para_ptr;
        abs_strength = vcl_abs(*strength_ptr);
        
        *orientation_ptr = *theta_ptr + half_pi;
        if (*orientation_ptr > two_pi)
          *orientation_ptr -= two_pi;
      }

      if (vcl_abs(*perp_ptr) > abs_strength)
      {
        *strength_ptr = *perp_ptr;
        *orientation_ptr = *theta_ptr;
        if (*orientation_ptr > two_pi)
          *orientation_ptr -= two_pi;
      }

      para_ptr += para_istep;
      perp_ptr += perp_istep;
      theta_ptr += theta_istep;
      strength_ptr += strength_istep;
      orientation_ptr += orientation_istep;
    }

    para_row_ptr += para_jstep;
    perp_row_ptr += perp_jstep;
    theta_row_ptr += theta_jstep;
    strength_row_ptr += strength_jstep;
    orientation_row_ptr += orientation_jstep;
  }
  
  // we could deal with scales outside this function
  // e.g. vil_gaussian_2nd_derivative_multiscale()
  // scale(swap_idx) = scales(ii);
}                          

#undef VIL_GAUSSIAN_DERIVATIVES_INSTANTIATE
#define VIL_GAUSSIAN_DERIVATIVES_INSTANTIATE(srcT, destT) \
template void vil_gaussian_1st_derivative(const vil_image_view<srcT >& src, \
                                 vil_image_view<destT >& line_strength, \
                                 vil_image_view<destT >& orientation, \
                                 double sd, int half_width, \
                                 vil_convolve_boundary_option bo); \
template void vil_gaussian_2nd_derivative(const vil_image_view<srcT >& src, \
                                 vil_image_view<destT >& line_strength, \
                                 vil_image_view<destT >& orientation, \
                                 double sd, int half_width, \
                                 vil_convolve_boundary_option bo)
                                 
#endif // __vil_gaussian_derivatives_txx__