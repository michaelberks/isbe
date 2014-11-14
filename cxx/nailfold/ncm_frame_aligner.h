#ifndef __ncm_frame_aligner_h__
#define __ncm_frame_aligner_h__

#include <vnl/vnl_matrix.h>

#include <vil/vil_image_view.h>

//: Class to align (register) two frames with respect to translation.
class ncm_frame_aligner
{
public: // Data types

  enum filter_type { filter_invalid = 0,
                     filter_g1d,
                     filter_g2d,
                     filter_sobel,
                     filter_first = filter_g1d,
                     filter_last = filter_sobel };

public: // methods  

  ncm_frame_aligner();

  //: Set the filtering method.
  void set_filter(filter_type filter);

  //: Return the filter name as a string.
  vcl_string filter_name();

  //: Set the number of levels of downsampling for efficiency.
  void set_n_levels(unsigned n_levels);

  //: Define pixels to ignore when matching 'interesting' points.
  void set_mask(vil_image_view<bool> mask_image);

  //: Return the current scaled mask image.
  vil_image_view<bool> mask() const;

  //: Set the image you want to move.
  void set_source(
    const vil_image_view<vxl_byte>& src);

  //: Return the source image.
  vil_image_view<vxl_byte> source_image();

  //: Set the image onto which you want to register src.
  void set_destination(
    const vil_image_view<vxl_byte>& dest);

  //: Return the destination image.
  vil_image_view<vxl_byte> destination_image();

  //: Return true if aligner has both a source and a destination image.
  bool is_ready();

  //: Compute the translation that maps the source onto the destination.
  void align_src_to_dest(  
    int di_init = 0,
    int dj_init = 0);

  //: Set the maximum displacement in x.
  void set_di_radius(
    int di_radius);

  //: Set the maximum displacement in y.
  void set_dj_radius(
    int dj_radius);

  //: Swap the 'pointers' to src and dest points.
  void swap_src_with_dest();

  //:
  void set_size(unsigned ni, unsigned nj);

  //: Size of the aligned images.
  unsigned ni();
  unsigned nj();

  //: Copy the likely displacements into two reference parameters.
  void get_displacements(
    int& peak_di, 
    int& peak_dj);

  //: Compute the linear photometric transformation needed to match src to dest
  //  i.e. (scale * src) + offset = dest
  void align_greylevels(
    double& scale, 
    double& offset);

  //: Return the mean squared error between overlapping parts of the images.
  double mean_squared_error();

private: // methods

  //:
  void get_interest_points(
    const vil_image_view<vxl_byte> image,
    vil_image_view<bool>& points_image);

  //:
  void get_interest_points_g2d(
    const vil_image_view<vxl_byte> image,
    vil_image_view<bool>& points_image);

  //:
  void get_interest_points_g1d(
    const vil_image_view<vxl_byte> image,
    vil_image_view<bool>& points_image);

  //:
  void get_interest_points_sobel(
    const vil_image_view<vxl_byte> image,
    vil_image_view<bool>& points_image);

  //:
  void count_votes(
    const vil_image_view<bool>& image,
    unsigned i, unsigned j, 
    vnl_matrix<unsigned>& accumulator);

  //:
  void find_accumulator_peak(
    const vnl_matrix<unsigned>& acc,
    int& peak_di, int& peak_dj);

  //:
  void output_to_log();

  //: Return the part of the source image that overlaps with the destination.
  vil_image_view<vxl_byte> src_intersection();

  //: Return the part of the destination image that overlaps with the source.
  vil_image_view<vxl_byte> dest_intersection();

  //: Reduce the mask image to correspond with the image pyramid level.
  void decimate_mask();

private: // constants

private: // variables

  //: Number of levels of subsampling in x and y (factor of two every time).
  unsigned n_levels_;

  //: Filtering method.
  filter_type filter_;

  //: Maximum displacement in x.
  int di_radius_;

  //: Maximum displacement in y.
  int dj_radius_;

  //: Position of the most likely offset
  int peak_di_;
  int peak_dj_;

  bool displacements_are_valid_;
  
  //: Source and destination images. The source image is translated by
  //  (peak_di_, peak_dj_) such that it is aligned with the destination.
  vil_image_view<vxl_byte> src_image_;
  vil_image_view<vxl_byte> dest_image_;

  //: Whether or not to apply a mask.
  bool use_mask_;

  //: Image that defines which pixels to ignore when matching features.
  vil_image_view<bool> mask_;
  vil_image_view<bool> mask_scaled_;

  //: Interesting points in source and destination images.
  //  Stored so that they can be reused for subsequent frames without 
  //  recomputing.
  vil_image_view<bool> src_pts_;
  vil_image_view<bool> dest_pts_;
};

#endif __ncm_frame_aligner_h__