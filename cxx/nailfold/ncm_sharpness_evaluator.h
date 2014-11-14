#ifndef __ncm_sharpness_evaluator__
#define __ncm_sharpness_evaluator__

#include <vil/vil_image_view.h>

class ncm_sharpness_evaluator
{

//  INTERFACE

public:

  //: Default ctor
  ncm_sharpness_evaluator();

  //: Return the sharpness of the image, as defined by the chosen method.
  double sharpness(const vil_image_view<vxl_byte>& vxl_image);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private: // Methods

  //: Define the region of interest of the image
  vil_image_view<vxl_byte> image_roi(
    const vil_image_view<vxl_byte>& vxl_image);

  //: Criterion based on Sobel filtering.
  double criterion_sobel(
    const vil_image_view<vxl_byte>& vxl_image);

  double criterion_laplacian(
    const vil_image_view<vxl_byte>& vxl_image);
  
  double criterion_entropy(
    const vil_image_view<vxl_byte>& vxl_image);
  
  double criterion_fft(
    const vil_image_view<vxl_byte>& vxl_image);

  double criterion_local_variance(
    const vil_image_view<vxl_byte>& vxl_image);

private: // Variables

  unsigned buffer_size_;
};

#endif // __ncm_sharpness_evaluator__