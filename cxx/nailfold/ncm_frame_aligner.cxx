#include "ncm_frame_aligner.h"

#include <vcl_vector.h>
#include <vcl_cmath.h>

#include <vil/vil_math.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_decimate.h>

#include <vil/algo/vil_gauss_reduce.h>
#include <vil/algo/vil_convolve_1d.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_suppress_non_max.h>
#include <vil/algo/vil_suppress_non_max_edges.h>
#include <vil/algo/vil_threshold.h>

#include <vnl/vnl_matrix.h>

#include <nailfold/vil_gaussian_derivatives.h>
#include <nailfold/vil_suppress_non_max_dir.h>

//
// Public member functions
//
 
//: Constructor
ncm_frame_aligner::ncm_frame_aligner()
: n_levels_(0),
  filter_(filter_sobel),
  di_radius_(40),
  dj_radius_(40),
  displacements_are_valid_(false),
  peak_di_(vcl_numeric_limits<int>::quiet_NaN()),
  peak_dj_(vcl_numeric_limits<int>::quiet_NaN()),
  mask_(vil_image_view<bool>(0,0)),
  use_mask_(false)
{
}
 
//: Set the filtering method.
void ncm_frame_aligner::set_filter(filter_type filter)
{
  if ((filter_first <= filter) && (filter <= filter_last))
    filter_ = filter;
}
 
//: Return the filter name as a string.
vcl_string ncm_frame_aligner::filter_name()
{
  switch (filter_)
  {
    case filter_g1d:
      return "g1d"; break;
    case filter_g2d:
      return "g2d"; break;
    case filter_sobel:
      return "sobel"; break;
  }

  return "<unknown>";
}
 
//: Set the number of levels of downsampling for efficiency.
void ncm_frame_aligner::set_n_levels(unsigned n_levels)
{
  n_levels_ = n_levels;
  if (use_mask_)
    decimate_mask();
}
 
//: Define pixels to ignore when matching 'interesting' points.
void ncm_frame_aligner::set_mask(vil_image_view<bool> mask_image)
{
  use_mask_ = true;
  mask_ = mask_image;
  decimate_mask();
}
 
//: Return the current scaled mask image.
vil_image_view<bool> ncm_frame_aligner::mask() const
{
  return mask_scaled_;
  //return mask_;
}
   
//: Set the image you want to move.
void ncm_frame_aligner::set_source(const vil_image_view<vxl_byte>& img)
{
  vil_image_view<vxl_byte> img1 = img;
  vil_image_view<vxl_byte> img2;

  src_image_ = img1;

  // Subsample the image n_levels_ times.
  if (n_levels_ > 0)
  {
    for (unsigned i = 0; i < n_levels_; ++i)
    {
      vil_gauss_reduce_121(img1, img2);
      vil_image_view<vxl_byte> tmp = img1;
      img1 = img2; // img1 = small
      img2 = tmp;  // img2 = big
    }
  }
  
  get_interest_points(img1, src_pts_);

  displacements_are_valid_ = false;
}
 
//: Return the source image.
vil_image_view<vxl_byte> ncm_frame_aligner::source_image()
{
  return src_image_;
}
 
//: Set the image onto which you want to register src.
void ncm_frame_aligner::set_destination(const vil_image_view<vxl_byte>& img)
{
  vil_image_view<vxl_byte> img1 = img;
  vil_image_view<vxl_byte> img2;

  dest_image_ = img1;

  // Subsample the image n_levels_ times.
  if (n_levels_ > 0)
  {
    for (unsigned i = 0; i < n_levels_; ++i)
    {
      vil_gauss_reduce_121(img1, img2);
      vil_image_view<vxl_byte> tmp = img1;
      img1 = img2; // img1 = small
      img2 = tmp;  // img2 = big
    }
  }

  get_interest_points(img1, dest_pts_);

  displacements_are_valid_ = false;
}
 
//: Return the destination image.
vil_image_view<vxl_byte> ncm_frame_aligner::destination_image()
{
  return dest_image_;
}
 
//: Return true if aligner has both a source and a destination image.
bool ncm_frame_aligner::is_ready()
{
  return ( (src_pts_.ni() > 0) && 
           (dest_pts_.ni() > 0) );
}
 
//: Return translation that maps src to dest.
void ncm_frame_aligner::align_src_to_dest(
  int di_init /* = 0 */,
  int dj_init /* = 0 */)
{

  // Downscale initial offset
  di_init = di_init >> n_levels_;
  dj_init = dj_init >> n_levels_;

  // Downscale displacement limits.
  const int di_radius = di_radius_ >> n_levels_;
  const int dj_radius = dj_radius_ >> n_levels_;

  // Add votes to accumulator for every feature pixel in source image.
  vnl_matrix<unsigned> accumulator(2*di_radius + 1,
                                   2*dj_radius + 1,
                                   0);

  const unsigned imin = di_radius + vcl_abs(di_init);
  const unsigned imax = src_pts_.ni() - imin;
  const unsigned jmin = dj_radius + vcl_abs(dj_init);
  const unsigned jmax = src_pts_.nj() - jmin;

  for (unsigned i = imin; i < imax; ++i)
  {
    for (unsigned j = jmin; j < jmax; ++j)
    {
      if (src_pts_(i,j) && (!use_mask_ || !mask_(i,j)))
        count_votes(dest_pts_, i+di_init, j+dj_init, accumulator);
    }
  }

  find_accumulator_peak(accumulator,
                        peak_di_, peak_dj_);
  peak_di_ += di_init;
  peak_dj_ += dj_init;

  // Upscale the estimated displacement.
  if (n_levels_ > 0)
  {
    peak_di_ = (peak_di_ << n_levels_);
    peak_dj_ = (peak_dj_ << n_levels_);
  }

  displacements_are_valid_ = true;
}
 
//: Set the maximum displacement in x.
void ncm_frame_aligner::set_di_radius(int di_radius)
{
  di_radius_ = vcl_abs(di_radius);
}
 
//: Set the maximum displacement in y.
void ncm_frame_aligner::set_dj_radius(int dj_radius)
{
  dj_radius_ = vcl_abs(dj_radius);
}
 
//: Swap the 'pointers' to src and dest points.
void ncm_frame_aligner::swap_src_with_dest()
{
  vil_image_view<bool> temp_pts_ = src_pts_;
  src_pts_ = dest_pts_;
  dest_pts_ = temp_pts_;

  vil_image_view<vxl_byte> temp_image_ = src_image_;
  src_image_ = dest_image_;
  dest_image_ = temp_image_;

  displacements_are_valid_ = false;
}
 
//:
void ncm_frame_aligner::set_size(unsigned ni, unsigned nj)
{
} 
 
//: Return the width of the incoming aligned images.
unsigned ncm_frame_aligner::ni()
{
  return src_image_.ni();
}
 
//: Return the height of the incoming aligned images.
unsigned ncm_frame_aligner::nj()
{
  return src_image_.nj();
}
 
//: 
void ncm_frame_aligner::get_displacements(
  int& peak_di, 
  int& peak_dj)
{
  if (displacements_are_valid_)
  {
    peak_di = peak_di_;
    peak_dj = peak_dj_;
  }
  else
  {
    peak_di_ = vcl_numeric_limits<int>::quiet_NaN();
    peak_dj_ = vcl_numeric_limits<int>::quiet_NaN();
  }
}
 
//: Compute the linear photometric transformation needed to match src to dest
//  i.e. (scale * src) + offset = dest
void ncm_frame_aligner::align_greylevels(
  double& scale, 
  double& offset)
{
  vil_image_view<double> src;
  vil_convert_cast(src_intersection(), src);
  double src_mean, src_var;
  vil_math_mean_and_variance(src_mean, src_var, src, 0);

  vil_image_view<double> dest;
  vil_convert_cast(dest_intersection(), dest);
  double dest_mean, dest_var;
  vil_math_mean_and_variance(dest_mean, dest_var, dest, 0);

  // Estimating the gain while using the image intersections drives the value
  // toward zero, washing out all contrast. Leave it out for now until I've
  // figured out why. (Standard deviation is heavily dependent on the ROI - much
  // more so than the mean.)
  scale = 1.0;
  //scale = vcl_sqrt(dest_var/src_var);

  offset = 0.0;
  offset = dest_mean - scale*src_mean;
}
 
//: Return the MSE in grey levels between the overlapping portions of the
//  aligned images.
double ncm_frame_aligner::mean_squared_error()
{
  vil_image_view<double> src;
  vil_convert_cast(src_intersection(), src);
  vil_math_normalise(src);

  vil_image_view<double> dest;
  vil_convert_cast(dest_intersection(), dest);
  vil_math_normalise(dest);

  const double type_indicator = 0.0;
  const double ssd = vil_math_ssd(src, dest, type_indicator);
  const unsigned n_pixels = src_intersection().ni() * src_intersection().nj();

  return ssd / n_pixels;
}

//
// Private member functions
//
 
//: Process an image to find a discrete set of interesting points.
void ncm_frame_aligner::get_interest_points(
  const vil_image_view<vxl_byte> image,
  vil_image_view<bool>& points_image)
{
  switch (filter_)
  {
  case filter_g1d:
    get_interest_points_g1d(image, points_image);
    break;

  case filter_g2d:
    get_interest_points_g2d(image, points_image);
    break;

  case filter_sobel:
  default:
    get_interest_points_sobel(image, points_image);
  }
}
 
//: Process an image to find a discrete set of interesting points from gaussian
//  second derivatives.
void ncm_frame_aligner::get_interest_points_g2d(
  const vil_image_view<vxl_byte> image,
  vil_image_view<bool>& points_image)
{
  const vil_convolve_boundary_option bo = vil_convolve_constant_extend;

  // Run the filters (first derivatives, second derivatives, both) over the 
  // images.
  vil_image_view<double> strength;
  vil_image_view<double> orientation;
  vil_gaussian_2nd_derivative(image,
                              strength, orientation,
                              4.0, -1, bo);

  // Process the filtered images (e.g. nonmaximal suppression) to get feature 
  // points.
  vil_image_view<double> maxima;
  vil_suppress_non_max_dir(strength, orientation,
                           maxima,
                           /* normals = */ false);

  // Suppress non-maximal points
  vil_math_rms(maxima, maxima);
  vil_threshold_above(maxima, points_image, 0.5);
}
 
//: Process an image to find a discrete set of interesting points from gaussian
//  first derivatives.
void ncm_frame_aligner::get_interest_points_g1d(
  const vil_image_view<vxl_byte> image,
  vil_image_view<bool>& points_image)
{
  const vil_convolve_boundary_option bo = vil_convolve_constant_extend;

  // Run the filters (first derivatives, second derivatives, both) over the 
  // images.
  vil_image_view<double> strength;
  vil_image_view<double> orientation;
  vil_gaussian_1st_derivative(image,
                              strength, orientation,
                              4.0, -1, bo);

  // Process the filtered images (e.g. nonmaximal suppression) to get feature 
  // points.
  vil_image_view<double> maxima;
  vil_suppress_non_max_dir(strength, orientation,
                           maxima,
                           /* normals = */ true);

  // Suppress non-maximal points
  vil_math_rms(maxima, maxima);
  vil_threshold_above(maxima, points_image, 0.5);
}
 
//: Process an image to find a discrete set of interesting points from Sobel
//  filtering.
void ncm_frame_aligner::get_interest_points_sobel(
  const vil_image_view<vxl_byte> image,
  vil_image_view<bool>& points_image)
{
  // Get gradients in x and y.
  vil_image_view<double> grad_ij;
  vil_sobel_3x3(image, grad_ij);

  vil_image_view<double> grad_magnitude;
  vil_math_rss(grad_ij, grad_magnitude);

  const double grad_threshold = 1.0;
  vil_threshold_above(grad_magnitude, points_image, grad_threshold);
}
 
//: Update the accumulator for displacements, given a putative centre.
void ncm_frame_aligner::count_votes(
  const vil_image_view<bool>& image,
  unsigned i, unsigned j, 
  vnl_matrix<unsigned>& accumulator)
{
  // Downscale displacement limits.
  const int di_radius = (di_radius_ >> n_levels_);
  const int dj_radius = (dj_radius_ >> n_levels_);

  const vcl_ptrdiff_t istep = image.istep();
  const vcl_ptrdiff_t jstep = image.jstep();
  
  int acc_r = 0;
  bool const* row_ptr = image.top_left_ptr() + 
                        (i-di_radius)*istep + 
                        (j-dj_radius)*jstep;

  for (int row = -dj_radius; row <= dj_radius; ++row)
  {
    unsigned* acc_ptr = accumulator[acc_r];

    bool const* col_ptr = row_ptr;
    for (int col = -di_radius; col <= dj_radius; ++col)
    {

      if (*col_ptr)
        ++(*acc_ptr);

      col_ptr += istep;
      ++acc_ptr;
    }

    row_ptr += jstep;
    ++acc_r;
  }
}
 
//: Find the displacement with maximal value in an accumulator.
void ncm_frame_aligner::find_accumulator_peak(
  const vnl_matrix<unsigned>& acc,
  int& peak_di,
  int& peak_dj)
{
  // Downscale displacement limits.
  const int di_radius = (di_radius_ >> n_levels_);
  const int dj_radius = (dj_radius_ >> n_levels_);

  unsigned max_votes = 0;
  
  peak_di = peak_dj = 0;

  for (int dj = -dj_radius; dj <= dj_radius; ++dj)
  {
    for (int di = -di_radius; di <= di_radius; ++di)
    {
      if (max_votes < acc(dj+dj_radius, di+di_radius))
      {
        max_votes = acc(dj+dj_radius, di+di_radius);
        peak_di = di;
        peak_dj = dj;
      }
    }
  }
}
 
//: 
void ncm_frame_aligner::output_to_log()
{
  //vcl_cout << "Writing output..." << vcl_endl;

  //// Output results.
  //vil_image_view<vxl_byte> imgout;
  //vcl_string outname;

  //vil_convert_cast(mosaic, imgout);
  //outname = rootdir + "_mosaic.png";
  //vil_save(imgout, outname.c_str());

  //vcl_ofstream ofs;
  //vcl_string displacement_log = rootdir + "_displacements.txt";
  //ofs.open(displacement_log.c_str());
  //{
  //  ofs << "Image: " << vcl_setw(4) << 1 
  //      << "  di_abs: " << vcl_setw(4) << di[0]
  //      << "  dj_abs: " << vcl_setw(4) << dj[0]
  //      << "  di_rel: " << vcl_setw(3) << 0 
  //      << "  dj_rel: " << vcl_setw(3) << 0 << vcl_endl;
  //  for (unsigned i = 1; i < n_images; ++i)
  //  {
  //    ofs << "Image: " << vcl_setw(4) << i+1
  //        << "  di_abs: " << vcl_setw(4) << di[i]
  //        << "  dj_abs: " << vcl_setw(4) << dj[i]
  //        << "  di_rel: " << vcl_setw(3) << di[i]-di[i-1]
  //        << "  dj_rel: " << vcl_setw(3) << dj[i]-dj[i-1] << vcl_endl;
  //  }
  //}
  //ofs.close();
}
 
//: Return the part of the source image that overlaps with the destination.
vil_image_view<vxl_byte> ncm_frame_aligner::src_intersection()
{
  if (displacements_are_valid_)
  {
    const int i0 = vcl_max(0, -peak_di_);
    const int ni = src_image_.ni() - vcl_abs(peak_di_);
    const int j0 = vcl_max(0, -peak_dj_);
    const int nj = src_image_.nj() - vcl_abs(peak_dj_);

    return vil_crop(src_image_, i0, ni, j0, nj);
  }
  else
    return vil_image_view<vxl_byte>(0,0);
}
 
//: Return the part of the destination image that overlaps with the source.
vil_image_view<vxl_byte> ncm_frame_aligner::dest_intersection()
{
  if (displacements_are_valid_)
  {
    const int i0 = vcl_max(0, peak_di_);
    const int ni = dest_image_.ni() - vcl_abs(peak_di_);
    const int j0 = vcl_max(0, peak_dj_);
    const int nj = dest_image_.nj() - vcl_abs(peak_dj_);

    return vil_crop(dest_image_, i0, ni, j0, nj);
  }
  else
    return vil_image_view<vxl_byte>(0,0);
}
 
//: Reduce the mask image to correspond with the image pyramid level.
void ncm_frame_aligner::decimate_mask()
{
  if ((mask_.ni() > 0) && (mask_.nj() > 0))
  {
    mask_scaled_ = mask_;

    for (unsigned L = 0; L < n_levels_; ++L)
    {
      vil_image_view<bool> cropped_mask = 
        vil_crop(mask_scaled_, 1, mask_scaled_.ni()-1,
                               1, mask_scaled_.nj()-1);
      mask_scaled_ = vil_decimate(cropped_mask, 2);
    }
  }
}