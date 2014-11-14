// vcl
#include <vcl_vector.h>
#include <vcl_string.h>
#include <vcl_algorithm.h>

// vxl/core
#include <vil/vil_image_view.h>
#include <vil/vil_load.h>

#include <vul/vul_arg.h>
#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>

#include <nailfold/flow/ncm_flow_field.h>
#include <nailfold/flow/ncm_flow_estimator_cholmod.h>

// vxl/contrib

// isbe
//#include <nailfold/ncm_flow/ncm_flow_brightness_cost.h>

enum error_code {
  SUCCESS = 0,
  INVALID_INPUT_PATH,
  INVALID_OUTPUT_PATH
};

error_code check_input_path(vcl_string& input_path)
{
  // Parse string to ensure it is in a sensible format

  if ( !vul_file::is_directory(input_path) )
    return INVALID_INPUT_PATH;

  return SUCCESS;
}

error_code check_output_path(vcl_string& output_path)
{
  // Parse string to ensure it is in a sensible format
  // Check that the string is a valid path
  // Create the path if necessary
  if (!vul_file::is_directory(output_path))
  {
    vul_file::make_directory_path(output_path);
  }

  return SUCCESS;
}

vcl_vector<vcl_string> files_in(
  vcl_string input_path)
{
  vcl_string search_str = input_path + "/frame*.png";

  vcl_vector<vcl_string> filenames;

  vul_file_iterator fn(search_str.c_str());
  for (fn; fn; ++fn)
    filenames.push_back(fn());

  // Ensure filenames are sorted (vul_file_iterator doesn't)
  vcl_sort(filenames.begin(), filenames.end());

  return filenames;
}

void load_image_sequence_from(
  const vcl_string& input_path,
  vcl_vector< vil_image_view<vxl_byte> >& image_stack,
  int n_images = -1)
{
  // Get list of images
  vcl_vector<vcl_string> filenames = files_in(input_path);

  if ( (0 <= n_images) && (n_images <= filenames.size()) )
    filenames.resize(n_images);

  n_images = filenames.size();
  image_stack.resize(n_images);

  // Read images in one by one
  for (unsigned i = 0; i < n_images; ++i)
  {
    image_stack[i] = vil_load(filenames[i].c_str());
  }
}
 
//: Main function
int main( int argc, char* argv[] )
{
  error_code ec;

  vul_arg<vcl_string> arg_input_path("-i",
                                     "Input image path", ".");
  vul_arg<unsigned> arg_n_images("-nimages",
                                 "Max. number of images", 60);
  vul_arg<unsigned> arg_n_levels("-nlevels",
                                 "Max. number of pyramid levels", 2);
  vul_arg_parse(argc, argv);

  vcl_string input_path = arg_input_path();
  ec = check_input_path(input_path);
  if (ec != SUCCESS)
    return ec;

  vcl_string output_path = input_path + "/flow.cpp/";
  ec = check_output_path(output_path);
  if (ec != SUCCESS)
    return ec;

  vcl_vector< vil_image_view<vxl_byte> > image_stack;
  load_image_sequence_from(input_path, image_stack, arg_n_images());
  
  // Minimize the error function
  ncm_flow_estimator_cholmod flow_estimator;
  flow_estimator.set_output_root(output_path);

  ncm_flow_field flow_field;
  flow_estimator.estimate_flow_from(image_stack, flow_field, arg_n_levels());

  return SUCCESS;
}