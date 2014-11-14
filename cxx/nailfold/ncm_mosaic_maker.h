#ifndef __ncm_mosaic_maker_h__
#define __ncm_mosaic_maker_h__

#include <vcl_vector.h>
#include <vcl_string.h>

#include <vil/vil_image_view.h>

class ncm_mosaic_maker
{
public:

  //: Constructor.
  ncm_mosaic_maker(
    unsigned input_ni,
    unsigned input_nj);

  //: Make a mosaic from a list of filenames and their associated relative
  //  displacements.
  bool make_mosaic_from(
    const vcl_vector<vcl_string>& filenames,
    const vcl_vector<int>& di,
    const vcl_vector<int>& dj);

	//: Set the maximum size (in number of pixels) for the mosaic - this stops the
	//mosaic from crashing the computer if we go crazy with the motors!
	void set_max_pixels(unsigned max);

  //: Set the size of the mosaic and its origin (w.r.t. the top-left pixel).
  void set_mosaic_limits(
    unsigned mosaic_ni,
    unsigned mosaic_nj,
    unsigned origin_i = 0,
    unsigned origin_j = 0);

  //: Add an image to the mosaic at a displacement of (di,dj) from the current
  //  location.
  virtual 
  void add_image(
    const vil_image_view<vxl_byte> img,
    unsigned di, unsigned dj);

  //; Return the final mosaic.
  vil_image_view<vxl_byte> mosaic();

  //: Return the map of frame counts. i.e. The value at a pixel indicates 
  //  the number of frames in which the pixel was visible.
  vil_image_view<unsigned> count_image();

protected: //Methods

	 //: Determine the size of the mosaic from a sequence of displacements.
  bool get_mosaic_limits(
    const vcl_vector<int>& di,
    const vcl_vector<int>& dj);

protected: // constants

  //: Size of the input images.
  const unsigned input_ni_;
  const unsigned input_nj_;

protected: // variables

  //: Location of the origin (i.e. image #1) in the mosaic.
  unsigned origin_i_;
  unsigned origin_j_;

private: // methods

  //: Private default constructor. Must give an input image size.
  ncm_mosaic_maker();

  void initialize_workspace();
   
  //: Create the map that weights every pixel for a smooth blend.
  void create_weight_image();
   
  //:
  void update_weights();
   
  //:
  void update_count();

  //:
  void add_to_mosaic(
    const vil_image_view<vxl_byte>& image);

private: // variables

  bool workspace_is_initialized_;
  bool mosaic_is_valid_;
    
  unsigned mosaic_ni_;
  unsigned mosaic_nj_;

  unsigned current_i_;
  unsigned current_j_;

	unsigned max_pixels_;

  //: Map of pixel weights. Equal in size to the input image.
  vil_image_view<double> weights_image_;

  //: The output mosaic.
  vil_image_view<double> mosaic_;

  //: Total weight at every pixel in the mosaic.
  vil_image_view<double> weights_;

  //: Total count at every pixel in the mosaic.
  vil_image_view<unsigned> count_;

  vil_image_view<double> mosaic_frame_;
  vil_image_view<double> weights_frame_;

  //: Input image, weighted by the map stored in weights_image_.
  vil_image_view<double> weighted_image_;

  //: The final mosaic in vxl_byte format.
  vil_image_view<vxl_byte> mosaic_byte_;
};

#endif // __ncm_mosaic_maker_h__