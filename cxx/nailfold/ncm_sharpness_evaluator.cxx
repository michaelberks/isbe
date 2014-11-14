#include "ncm_sharpness_evaluator.h"

#include <vil/vil_plane.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>
#include <vil/algo/vil_threshold.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_convolve_2d.h>
#include <vil/algo/vil_histogram.h>
#include <vil/algo/vil_fft.h>

//
// Image sharpness (focus) criterion functions.
// See Krotkov, "Focusing", IJCV (1), 1987.
//

//
//: Default constructor
ncm_sharpness_evaluator::ncm_sharpness_evaluator()
: buffer_size_(1)
{
}

//
//: Proxy for sharpness function, including a running average.
double ncm_sharpness_evaluator::sharpness(
    const vil_image_view<vxl_byte>& vxl_image)
{
  vil_image_view<vxl_byte> roi = image_roi(vxl_image);
  if (roi.ni() == 0)
    return -1.0;

  // Use static variables here because they are only relevant to this function 
  // and no other, so no need to make them members of the class.
  static vcl_vector<double> sharpness_values(buffer_size_, -1.0);
  static unsigned sharpness_index = 0;

  //sharpness_values[sharpness_index] = criterion_sobel(roi);
  sharpness_values[sharpness_index] = criterion_local_variance(roi);

  // Update ring buffer index.
  ++sharpness_index;
  if (sharpness_index == buffer_size_)
    sharpness_index = 0;

  // Compute running average.
  double mean_sharpness = 0.0;
  for (unsigned i = 0; i < buffer_size_; ++i)
    mean_sharpness += sharpness_values[i];
  mean_sharpness /= buffer_size_;

  return mean_sharpness;
}

//
// Private methods
//

//
//: Get a region of interest from the image.
//  Could be used to evaluate sharpness on only a portion of the image for
//  efficiency.
vil_image_view<vxl_byte> ncm_sharpness_evaluator::image_roi(
    const vil_image_view<vxl_byte>& vxl_image)
{
  // Set width and height of the region of interest
  unsigned ni = vxl_image.ni() / 2;
  unsigned nj = vxl_image.nj() / 2;

  unsigned i0 = (vxl_image.ni() - ni) / 2;
  unsigned j0 = (vxl_image.nj() - nj) / 2;

  return vil_crop(vil_plane(vxl_image, 0), i0, ni, j0, nj);
}
 
//
//: Sobel-based criterion
//  Works pretty well.
double ncm_sharpness_evaluator::criterion_sobel(
    const vil_image_view<vxl_byte>& vxl_image)
{
  // Compute x and y gradients
  vil_image_view<float> sobel;
  vil_sobel_3x3(vxl_image, sobel);
  
  // Compute edge strength
  vil_image_view<float> edge_strength;
  vil_math_image_vector_mag(vil_plane(sobel, 0), vil_plane(sobel, 1),
                            edge_strength);

  // Find strong edges
  vil_image_view<bool> strong_edges;
  vil_threshold_above(edge_strength, strong_edges, /* threshold = */ 0.0f);

  // Mask out weak edges
  vil_image_view<float> strong_edges_values;
  vil_convert_cast(strong_edges, strong_edges_values);

  vil_image_view<float> masked_edge_strength;
  vil_math_image_product(edge_strength, strong_edges_values, 
                         masked_edge_strength);

  // Sum over strong edges
  double edge_score = 0.0;
  //vil_math_sum(edge_score, masked_edge_strength, /* plane = */ 0);

  float min_val = 0, max_val = 0;
  vil_math_value_range(masked_edge_strength, min_val, max_val);
  edge_score = max_val;

  //double grey_sum = 0.0;
  //vil_math_sum(grey_sum, vxl_image, /* plane = */ 0);

  double grey_mean, grey_var;
  vil_math_mean_and_variance(grey_mean, grey_var, vxl_image, /* plane = */ 0);

  // Normalize with respect to overall grey level
  //return edge_sum / grey_sum;
  //return edge_sum / grey_var;
  return edge_score;
}

//
//: Laplacian-based criterion
//  Pants
double ncm_sharpness_evaluator::criterion_laplacian(
    const vil_image_view<vxl_byte>& vxl_image)
{
  // Compute approximation to Laplacian
  vil_image_view<float> laplacian;
  vil_image_view<float> kernel(3,3);
  float one_sixth = 1.0f / 6.0f;
  kernel(0,0) = kernel(0,2) = kernel(2,0) = kernel(2,2) = one_sixth;
  kernel(1,0) = kernel(0,1) = kernel(1,2) = kernel(2,1) = 4 * one_sixth;
  kernel(1,1) = -20 * one_sixth;

  float ac = 0.0;
  vil_convolve_2d(vxl_image, laplacian, kernel, ac);
  
  // Find strong responses
  vil_image_view<bool> strong_responses;
  vil_threshold_above(laplacian, strong_responses, /* threshold = */ 10.0f);

  // Mask out weak edges
  vil_image_view<float> strong_responses_values;
  vil_convert_cast(strong_responses, strong_responses_values);

  vil_image_view<float> masked_responses;
  vil_math_image_product(laplacian, strong_responses_values, 
                         masked_responses);

  // Sum over strong edges
  double criterion;
  vil_math_sum(criterion, masked_responses, /* plane = */ 0);

  return criterion;
}

//
//: Entropy-based criterion
//  Pants: only works if you can get really sharp images. Otherwise, the
//  distribution tends to the average value as blurring becomes extreme, giving
//  a sharply peaked distribution with low entropy.
double ncm_sharpness_evaluator::criterion_entropy(
    const vil_image_view<vxl_byte>& vxl_image)
{
  unsigned n_pixels = vxl_image.ni() * vxl_image.nj();
  vcl_vector<double> histo(256, 0);
  vil_histogram_byte(vxl_image, histo);

  double entropy = 0.0;
  for (unsigned i = 0; i < histo.size(); ++i)
  {
    if ((histo[i] != 0) &&      // p_i = 0
        (histo[i] != n_pixels)) // log(p_i) = 0
    {
      double p_i = histo[i] / n_pixels;
      entropy -= p_i * vcl_log(p_i);
    }
  }

  return -entropy; // minimize entropy rather than maximize it
}

//
//: FFT-based criterion
//  Theoretically, should be one of the best measures but determining what
//  exactly to extract from the FFT will have a big effect.
double ncm_sharpness_evaluator::criterion_fft(
    const vil_image_view<vxl_byte>& vxl_image)
{
  vil_image_view< vcl_complex<float> > fft_image;
  vil_convert_cast(vxl_image, fft_image);

  vil_fft_2d_fwd(fft_image);

  return 0.0;
}

//: Measure local variance in a region around every pixel and return some
//  property of the distribution (e.g. the mean or max).
double ncm_sharpness_evaluator::criterion_local_variance(
    const vil_image_view<vxl_byte>& vxl_image)
{
  vil_image_view<vxl_uint_32> int_image;
  vil_image_view<vxl_uint_32> int_squared_image;

  vil_math_integral_sqr_image(vxl_image,
                              int_image, int_squared_image);

  const int hw = 15; // half width

  const int fw = 2*hw + 1; // full_width
  const int np = fw*fw; // number of pixels

  vil_image_view<double> var_image(vxl_image.ni() - (2*hw),
                                   vxl_image.nj() - (2*hw));

  const unsigned ni = int_image.ni();
  const unsigned nj = int_image.nj();

  double max_var = -1.0;
  double sum_var = 0.0;
  unsigned count = 0;

  for (unsigned i = fw; i < ni; ++i)
  {
    for (unsigned j = fw; j < nj; ++j)
    {
      const long sum =   int_image(i, j)
                       - int_image(i-fw, j)
                       - int_image(i, j-fw)
                       + int_image(i-fw, j-fw);

      const double mean = static_cast<double>(sum) / np;

      const long sum_squared =   int_squared_image(i, j)
                               - int_squared_image(i-fw, j)
                               - int_squared_image(i, j-fw)
                               + int_squared_image(i-fw, j-fw);

      const double mean_squared = static_cast<double>(sum_squared) / np;
      const double var = mean_squared - (mean * mean);

      var_image(i-fw, j-fw) = var;

      if (var > max_var)
        max_var = var;

      if (var > 0.0)
      {
        ++count;
        sum_var += var;
      }
    }
  }

  const double local_var = sum_var / count;
  
  //double image_mean, image_var;
  //vil_math_mean_and_variance(image_mean, image_var, 
  //                           vxl_image, /* plane = */ 0);
  //return local_var / image_var;

  return local_var;
}
