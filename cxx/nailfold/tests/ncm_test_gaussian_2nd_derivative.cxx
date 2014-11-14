#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>

#include <nailfold/vil_gaussian_derivatives.h>

//
// Main function
int main( int argc, char* argv[] )
{
  // add your stuff here
  vcl_string filename = 
    "U:/projects/nailfold/images/oneloop.png";

  vil_image_view<vxl_byte> img = vil_load(filename.c_str());

  vil_image_view<float> line_strength;
  vil_image_view<float> orientation;

  // gaussian 2nd derivative processed image
  double std_dev = 3.0;
  unsigned half_width = static_cast<unsigned>(4*std_dev);
  vil_gaussian_2nd_derivative(img, 
                              line_strength, orientation, 
                              std_dev, half_width );

  vil_image_view<vxl_byte> b_image;
  vil_convert_stretch_range(line_strength, b_image);
  vil_save(b_image,"u:/projects/nailfold/images/ncm_test_strength.png");
  vil_convert_stretch_range(orientation, b_image);
  vil_save(b_image,"u:/projects/nailfold/images/ncm_test_orientation.png");

  return 0;
}


