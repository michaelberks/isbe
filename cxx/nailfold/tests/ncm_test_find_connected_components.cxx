#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/algo/vil_threshold.h>

//#include <nailfold/vil_find_connected_components.h>
#include <nailfold/vil_connected_components.h>

void threshold_sizes(const vil_image_view<int>& labels_image,
                     const vcl_vector<unsigned>& size_vector,
                     vil_image_view<bool>& bool_image)
{
  bool_image.set_size(labels_image.ni(), labels_image.nj());

  const unsigned size_threshold = 15;

  const int* lbl_row_ptr = labels_image.top_left_ptr();
  bool* out_row_ptr = bool_image.top_left_ptr();

  const unsigned ni = labels_image.ni();
  const unsigned nj = labels_image.nj();

  for (unsigned j = 0; j < nj; j++)
  {
    const int* lbl_ptr = lbl_row_ptr;
    bool* out_ptr = out_row_ptr;

    for (unsigned i = 0; i < ni; i++)
    {
      if (*lbl_ptr > 0)
        *out_ptr = (size_vector[(*lbl_ptr)-1] >= size_threshold);
      else
        *out_ptr = false;

      ++lbl_ptr;
      ++out_ptr;
    }

    lbl_row_ptr += labels_image.jstep();
    out_row_ptr += bool_image.jstep();
  }
}

//
// Main function
int main( int argc, char* argv[] )
{
  // add your stuff here
  vcl_string filename = 
    //"u:/projects/nailfold/images/peaks_image.png";
    "u:/projects/nailfold/images/conn_comp_test.png";

  vil_image_view<vxl_byte> img = vil_load(filename.c_str());
  vil_image_view<vxl_byte> flat_image;
  if (img.nplanes() > 1)
    vil_convert_planes_to_grey(img, flat_image);
  else
    flat_image.deep_copy(img);

  vil_image_view<bool> bool_image(flat_image.ni(), flat_image.nj());
  vil_threshold_above<vxl_byte>(flat_image, bool_image, 127);

  // Create connected component class based on bool_image
  vil_connected_components conn_comp(bool_image);

  // convert to vxl_byte and save
  vil_image_view<vxl_byte> b_image;
  vil_convert_stretch_range(conn_comp.label_image(), b_image);
  vil_save(b_image,"u:/projects/nailfold/images/conn_comp_output.png");

  // update binary image to include only those components with size greater
  // than a defined threshold
  vil_image_view<bool> size_image(flat_image.ni(), flat_image.nj());
  conn_comp.set_binning_strategy(vil_connected_components::BinningEqualCounts);
  vil_threshold_above<vxl_byte>(conn_comp.size_image(), size_image, 240);
  vil_convert_stretch_range(size_image, b_image);
  vil_save(b_image,"u:/projects/nailfold/images/conn_comp_thresholded.png");

  return 0;
}


