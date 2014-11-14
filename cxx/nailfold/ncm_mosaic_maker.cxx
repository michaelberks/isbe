#include "ncm_mosaic_maker.h"

#include <vil/vil_load.h>
#include <vil/vil_crop.h>
#include <vil/vil_math.h>
#include <vil/vil_convert.h>


//
// Public methods
//

ncm_mosaic_maker::ncm_mosaic_maker(
  unsigned input_ni,
  unsigned input_nj)
: // Constants
  input_ni_(input_ni),
  input_nj_(input_nj),
	max_pixels_(1e7),
  // Variables
  workspace_is_initialized_(false),
  mosaic_is_valid_(false),
  mosaic_ni_(0),
  mosaic_nj_(0),
  mosaic_(0,0),
  weights_(0,0),
  count_(0,0)
{
  create_weight_image();
}
 
//: Create a mosaic from a sequence of images and their offsets.
bool ncm_mosaic_maker::make_mosaic_from(
  const vcl_vector<vcl_string>& filenames,
  const vcl_vector<int>& di,
  const vcl_vector<int>& dj)
{
  if(!get_mosaic_limits(di, dj))
		return false;

  // Get displacements between successive images.
  const unsigned n_images = filenames.size();
  for (unsigned i = 0; i < n_images; ++i)
  {
    vil_image_view<vxl_byte> image = vil_load(filenames[i].c_str());
    add_image(image, di[i], dj[i]);
  }
	return true;
}
 
//: Define the maximum number of pixels the mosaic can contain
void ncm_mosaic_maker::set_max_pixels(unsigned max)
{
	max_pixels_ = max;
}

//: Define the size of the mosaic and its origin.
void ncm_mosaic_maker::set_mosaic_limits(
  unsigned mosaic_ni,
  unsigned mosaic_nj,
  unsigned origin_i /* = 0 */,
  unsigned origin_j /* = 0 */)
{
  mosaic_ni_ = mosaic_ni;
  mosaic_nj_ = mosaic_nj;
  origin_i_ = origin_i;
  origin_j_ = origin_j;

  current_i_ = origin_i_;
  current_j_ = origin_j_;

  initialize_workspace();
}

// virtual
void ncm_mosaic_maker::add_image(
  const vil_image_view<vxl_byte> img,
  unsigned di, unsigned dj)
{
  current_i_ += di;
  current_j_ += dj;

  // TODO: Normalize image with respect to existing mosaic grey values.

  add_to_mosaic(img);

  mosaic_is_valid_ = false;
}
 
//: Return the final result.
vil_image_view<vxl_byte> ncm_mosaic_maker::mosaic()
{
  if (!mosaic_is_valid_)
  {
    // Normalize with respect to total weights.
    vil_image_view<double> normalized_mosaic(mosaic_ni_, mosaic_nj_);
    vil_math_image_ratio(mosaic_, weights_, normalized_mosaic);

    vil_convert_cast(normalized_mosaic, mosaic_byte_);

    mosaic_is_valid_ = true;
  }

  return mosaic_byte_;
}
 
//: Return the map of frame counts. i.e. The value at a pixel indicates 
//  the number of frames in which the pixel was visible.
vil_image_view<unsigned> ncm_mosaic_maker::count_image()
{
  return count_;
}

//
// Private methods
//
 
void ncm_mosaic_maker::initialize_workspace()
{
  // Create mosaic image of desired size
  mosaic_.set_size(mosaic_ni_, mosaic_nj_);
  mosaic_.fill(0.0);

  weights_.set_size(mosaic_ni_, mosaic_nj_);
  weights_.fill(0.0);

  count_.set_size(mosaic_ni_, mosaic_nj_);
  count_.fill(0);

  workspace_is_initialized_ = true;
}
 
//: Determine the size of the mosaic from a sequence of displacements.
bool ncm_mosaic_maker::get_mosaic_limits(
  const vcl_vector<int>& displacement_x,
  const vcl_vector<int>& displacement_y)
{
  const unsigned n_images = displacement_x.size();

  assert(displacement_y.size() == n_images);

  if (n_images == 0)
    return false;

  // Compute cumulative displacements.
  // (i.e. displacements from a common datum.)
  vcl_vector<int> position_x(n_images, 0);
  vcl_vector<int> position_y(n_images, 0);
  for (unsigned i = 1; i < n_images; ++i)
  {
    position_x[i] = position_x[i-1] + displacement_x[i];
    position_y[i] = position_y[i-1] + displacement_y[i];
  }

  // Find minimum displacement in x and y.
  int min_x = 32767;
  int min_y = min_x;
  int max_x = -32767;
  int max_y = max_x;
  for (unsigned i = 0; i < n_images; ++i)
  {
    if (position_x[i] < min_x)
      min_x = position_x[i];
    if (max_x < position_x[i])
      max_x = position_x[i];

    if (position_y[i] < min_y)
      min_y = position_y[i];
    if (max_y < position_y[i])
      max_y = position_y[i];
  }

	unsigned mosaic_ni = (max_x - min_x) + input_ni_ + 1;
  unsigned mosaic_nj = (max_y - min_y) + input_nj_ + 1;

	//Check we've not tried to make to big a mosaic
	if ((mosaic_ni * mosaic_nj) > max_pixels_)
		return false;

	set_mosaic_limits(mosaic_ni, mosaic_nj,
										-min_x, -min_y);
	return true;

}
 
//: Create the map that weights every pixel for a smooth blend.
void ncm_mosaic_maker::create_weight_image()
{
  weights_image_.set_size(input_ni_, input_nj_);

  double mid_i = input_ni_/2;
  double mid_j = input_nj_/2;
  for (unsigned i = 0; i < input_ni_; ++i)
  {
    double weight_i = 1 - (vcl_abs(i - mid_i) / mid_i);

    for (unsigned j = 0; j < input_nj_; ++j)
    {
      double weight_j = 1 - (vcl_abs(j - mid_j) / mid_j);
      weights_image_(i,j) = weight_i * weight_j;
    }
  }
}
 
//: Update the map of pixel counts.
void ncm_mosaic_maker::update_count()
{
  vil_image_view<unsigned> view = vil_crop(count_, current_i_, input_ni_, 
                                                   current_j_, input_nj_);

  // Add one to the current view.
  vil_math_scale_and_offset_values(view, /* scale = */ 1, 
                                         /* offset = */ 1);
}
  
//: Update the map of total weight at every pixel.
void ncm_mosaic_maker::update_weights()
{
  vil_image_view<double> view = vil_crop(weights_, current_i_, input_ni_, 
                                                   current_j_, input_nj_);
  vil_math_image_sum(view, weights_image_, view);
}
 
//:
void ncm_mosaic_maker::add_to_mosaic(
  const vil_image_view<vxl_byte>& image)
{
  update_count();
  update_weights();

  // Weight the input image.
  vil_math_image_product(weights_image_, image, 
                         weighted_image_);

  // Add weighted image to mosaic.
  vil_image_view<double> view = vil_crop(mosaic_, current_i_, input_ni_, 
                                                  current_j_, input_nj_);
  vil_math_image_sum(view, weighted_image_, view);
}
