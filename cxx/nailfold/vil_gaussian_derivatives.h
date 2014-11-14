#ifndef __vil_gaussian_derivatives_h__
#define __vil_gaussian_derivatives_h__

#include <vil/vil_image_view.h>

#include <vil/algo/vil_convolve_1d.h>

template<class srcT, class destT>
void vil_gaussian_1st_derivative(const vil_image_view<srcT>& src,
                                 vil_image_view<destT>& line_strength,
                                 vil_image_view<destT>& orientation,
                                 double sd = 3.0, 
                                 int half_width = -1,
                                 vil_convolve_boundary_option bo = vil_convolve_zero_extend);
//void vil_gaussian_1st_derivative(const vil_image_view<vxl_byte>& src,
//                                 vil_image_view<double>& line_strength,
//                                 vil_image_view<double>& orientation,
//                                 double sd = 3.0, 
//                                 int half_width = -1);

template<class srcT, class destT>
void vil_gaussian_2nd_derivative(const vil_image_view<srcT>& src,
                                 vil_image_view<destT>& line_strength,
                                 vil_image_view<destT>& orientation,
                                 double sd = 3.0, 
                                 int half_width = -1,
                                 vil_convolve_boundary_option bo = vil_convolve_zero_extend);
//void vil_gaussian_2nd_derivative(const vil_image_view<vxl_byte>& src,
//                                 vil_image_view<double>& line_strength,
//                                 vil_image_view<double>& orientation,
//                                 double sd = 3.0, 
//                                 int half_width = -1);

#endif // __vil_gaussian_derivatives_h__