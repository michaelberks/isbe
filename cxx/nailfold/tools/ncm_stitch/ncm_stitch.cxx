#include <vcl_iomanip.h>
#include <vcl_algorithm.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_threshold.h>

#include <vil/algo/vil_gauss_reduce.h>

#include <vsl/vsl_quick_file.h>

#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>
#include <vul/vul_arg.h>

#include <nailfold/ncm_frame_aligner.h>
#include <nailfold/ncm_mosaic_maker.h>
#include <nailfold/ncm_sharpness_evaluator.h>
#include <nailfold/ncm_timebar.h>

//
//: Main function.
int main( int argc, char* argv[] )
{
  vul_arg<vcl_string> arg_root_dir("-i", "Input image path", ".");
  vul_arg<bool> arg_select("-select", "Select sharp images");
  vul_arg_parse(argc, argv);

  vcl_string rootdir = arg_root_dir();
  vcl_string search_str = rootdir + "/frame*.png";

  // Get list of files
  vul_file_iterator fn(search_str.c_str());
  vcl_vector<vcl_string> all_filenames;
  for (fn; fn; ++fn)
    all_filenames.push_back(fn());

  // Put files in alphabetical (and hopefully chronological) order.
  vcl_sort(all_filenames.begin(), all_filenames.end());

  const unsigned n = 400;
  if (all_filenames.size() > n)
    all_filenames.resize(n);

  unsigned n_images = all_filenames.size();
  if (n_images == 0)
    return 1;

  vil_image_view<vxl_byte> img1, img2;
  vcl_vector<double> all_sharpness_vec(n_images, -1.0);

  vil_image_view<double> mean_image;

  ncm_timebar tb("Computing sharpness and mask", n_images);
  ncm_sharpness_evaluator sharpness_evaluator;
  for (unsigned i = 0; i < n_images; ++i)
  {
    img1 = vil_load(all_filenames[i].c_str());
    all_sharpness_vec[i] = sharpness_evaluator.sharpness(img1);

    if (mean_image.ni() == 0)
      vil_convert_cast(img1, mean_image);
    else
    {
      vil_image_view<double> double_image;
      vil_convert_cast(img1, double_image);
      vil_math_image_sum(mean_image, double_image, mean_image);
    }

    tb.advance();
  }
  vil_math_scale_values(mean_image, 1.0/n_images);

  // Set the mask image
  vil_image_view<double> grad_ij;
  vil_sobel_3x3(mean_image, grad_ij);
  vil_image_view<double> grad_magnitude;
  vil_math_rss(grad_ij, grad_magnitude);
  double grad_threshold;
  vil_math_mean(grad_threshold, grad_magnitude, 0);
  vil_image_view<bool> mask_image;
  vil_threshold_above(grad_magnitude, mask_image, 2.0*grad_threshold);

  vcl_vector<vcl_string> filenames(0);
  vcl_vector<double> sharpness_vec(0);
  if (arg_select())
  {
    tb.reset("Selecting images", n_images-2);

    // Filter out those images that are less sharp than their neighbours.
    // Flag those for deletion by setting sharpness to -1.
    for (unsigned i = 1; i < n_images-1; ++i)
    {
      if ( (all_sharpness_vec[i] > all_sharpness_vec[i-1]) &&
           (all_sharpness_vec[i] > all_sharpness_vec[i+1]) )
      {
        sharpness_vec.push_back(all_sharpness_vec[i]);
        filenames.push_back(all_filenames[i]);
      }
      tb.advance();
    }
  }
  else
  {
    filenames = all_filenames;
    sharpness_vec = all_sharpness_vec;
  }

  n_images = filenames.size();

  unsigned input_ni = img1.ni();
  unsigned input_nj = img1.nj();


  tb.reset("Aligning images", n_images-1);

  ncm_frame_aligner aligner;

  const unsigned max_distance = 32;
  aligner.set_di_radius(max_distance);
  aligner.set_dj_radius(max_distance);

  aligner.set_n_levels(1);
  
  //aligner.set_filter(ncm_frame_aligner::filter_g1d);
  
  //aligner.set_mask(mask_image);

  //vil_image_view<vxl_byte> byte_image;
  //vil_convert_cast(aligner.mask(), byte_image);
  //vil_math_scale_values(byte_image, 255);
  //vil_save(byte_image, "u:/tmp/mask_image.png");

  // Load first valid image.
  img1 = vil_load(filenames[0].c_str());
  aligner.set_destination(img1);

  // Get displacements between successive images.
  vcl_vector<int> di(n_images, 0);
  vcl_vector<int> dj(n_images, 0);
  for (unsigned i = 1; i < n_images; ++i)
  {
    img1 = vil_load(filenames[i].c_str());
    aligner.set_source(img1);
    
    // Get displacement of img1 with respect to img2 (i.e. how img2 should
    // be translated to align the two).
    aligner.align_src_to_dest();
    aligner.get_displacements(di[i], dj[i]);

    aligner.swap_src_with_dest();

    tb.advance();
  }


  tb.reset("Creating mosaic", 1);
  ncm_mosaic_maker maker(input_ni, input_nj);
  maker.make_mosaic_from(filenames, di, dj);
  tb.advance();

  
  // Output results.
  tb.reset("Writing output", n_images);

  vcl_string output_dir = rootdir + "/mosaic/";
  vul_file::make_directory(output_dir);

  vcl_string select_suffix = "";
  if ( arg_select() )
    select_suffix = "_selected";

  vcl_string outname = output_dir + "mosaic_" + 
                       aligner.filter_name() + select_suffix + ".png";
  vil_save(maker.mosaic(), outname.c_str());

  vcl_ofstream ofs;
  vcl_string displacement_log = output_dir + "displacements.txt";
  ofs.open(displacement_log.c_str());
  {
    ofs << "Image: " << vcl_setw(4) << 1 
        << "  di_abs: " << vcl_setw(4) << di[0]
        << "  dj_abs: " << vcl_setw(4) << dj[0]
        << "  di_rel: " << vcl_setw(3) << 0 
        << "  dj_rel: " << vcl_setw(3) << 0 
        << "  sharpness: " << sharpness_vec[0]
        << vcl_endl;
    tb.advance();

    for (unsigned i = 1; i < n_images; ++i)
    {
      ofs << "Image: " << vcl_setw(4) << i+1
          << "  di_abs: " << vcl_setw(4) << di[i]+di[i-1]
          << "  dj_abs: " << vcl_setw(4) << dj[i]+dj[i-1]
          << "  di_rel: " << vcl_setw(3) << di[i]
          << "  dj_rel: " << vcl_setw(3) << dj[i] 
          << "  sharpness: " << sharpness_vec[i]
          << vcl_endl;
      tb.advance();
    }
  }
  ofs.close();

  return 0;
}
