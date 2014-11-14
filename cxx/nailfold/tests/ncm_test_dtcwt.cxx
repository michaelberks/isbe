#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_decimate.h>

#include <nailfold/vil_dtcwt.h>

template <class T>
void disp(const vil_image_view<T>& img)
{
  for (unsigned j = 0; j < img.nj(); ++j)
  {
    for (unsigned i = 0; i < img.ni(); ++i)
      vcl_cout << double(img(i,j)) << '\t';
    vcl_cout << vcl_endl;
  }
  vcl_cout << vcl_endl;
}

void disp_im(const vil_image_view< vcl_complex<double> >& img)
{
  for (unsigned j = 0; j < img.nj(); ++j)
  {
    for (unsigned i = 0; i < img.ni(); ++i)
      vcl_cout << img(i,j) << '\t';
    vcl_cout << vcl_endl;
  }
  vcl_cout << vcl_endl;
}

//
// Main function
int main( int argc, char* argv[] )
{
  unsigned N = 200;
  unsigned n_levels = 6;
  
  vil_image_view<vxl_byte> test_image(N,N);
  for (unsigned j = 0; j < test_image.nj(); ++j)
    for (unsigned i = 0; i < test_image.ni(); ++i)
      test_image(i,j) = ((i+2)*(j+1))%256;

  vil_dtcwt dt(test_image,n_levels);

  disp_im(vil_plane(dt.tree()[n_levels-1],0));

  return 0;
}


