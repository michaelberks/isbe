#include <vcl_cstddef.h>
#include <vcl_ctime.h>
#include <vcl_iosfwd.h>
#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_vector.h>
#include <vcl_complex.h>
#include <vcl_algorithm.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>

#include <vul/vul_file_iterator.h>

#include <nailfold/vil_dtcwt.h>

#include <nailfold/ncm_image_server.h>
#include <nailfold/ncm_scp_handler.h>
#include <nailfold/ncm_sharpness_evaluator.h>

//
// Main function
int main( int argc, char* argv[] )
{
  ncm_sharpness_evaluator sharpie;

  vcl_string root_dir = 
      "U:/projects/nailfold/capture/2013_03_06/Left/Digit4/x300/autofocus_vessels/";
  vcl_string search_str = "frame_*.png";

  // Get filenames.
  vul_file_iterator fn(root_dir + search_str);
  vcl_vector<vcl_string> filenames;
  for (fn; fn; ++fn)
  {
    filenames.push_back(fn());
  }
  vcl_sort(filenames.begin(), filenames.end());
  
  vcl_time_t t0 = vcl_clock();

  vcl_ofstream ofs;
  vcl_string logname = root_dir + "sharpness_test.log";
  ofs.open(logname.c_str());
  vcl_cout << "Evaluating images";
  for (unsigned i = 0; i < filenames.size(); ++i)
  {
    vil_image_view<vxl_byte> image = vil_load(filenames[i].c_str());
    ofs << sharpie.sharpness(image) << vcl_endl;
    vcl_cout << ".";
  }
  ofs.close();
  vcl_cout << vcl_endl;

  vcl_cout << "Completed in "
           << double(vcl_clock() - t0) / CLOCKS_PER_SEC 
           << " seconds." << vcl_endl;

  return 0;
}


