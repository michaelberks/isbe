#include "ncm_flow_estimator.h"

#include <vcl_sstream.h>
#include <vcl_ctime.h>

#include <vil/vil_save.h>
#include <vil/vil_pyramid_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_bilin_interp.h>

#include <vul/vul_file.h>

#include <nailfold/flow/ncm_flow_field.h>
#include <nailfold/flow/ncm_flow_null_robustifier.h>

//
//  Public methods
//
 
//:
ncm_flow_estimator::ncm_flow_estimator()
: output_root_("./")
{
}

void ncm_flow_estimator::set_output_root(
  vcl_string output_root)
{
  // Parse string to ensure it is in a sensible format

  // Check that the string is a valid path

  // Create the path if necessary
  if (!vul_file::is_directory(output_root))
    vul_file::make_directory_path(output_root);

  output_root_ = output_root;
}
 
//: Estimate flow from an image stack with n_levels levels.
void ncm_flow_estimator::estimate_flow_from(
  // Inputs
  vcl_vector< vil_image_view<vxl_byte> > const& image_stack,
  // Outputs
  ncm_flow_field& flow_field,
  // Optional inputs
  unsigned n_levels /* = 1 */)
{
  vcl_vector< vcl_vector< vil_image_view<float> > > image_stack_pyramid;
  create_image_pyramid_from(image_stack, n_levels, 
                            image_stack_pyramid);

  flow_field.set_size(image_stack_pyramid[n_levels-1][1].ni(),
                      image_stack_pyramid[n_levels-1][1].nj());

  vcl_vector< vil_image_view<float> > warped_image_stack;

  for (int L = n_levels-1; L >= 0; --L)
  {
    time_t start_time = time(NULL);
    estimate_flow_at_level(image_stack_pyramid[L], warped_image_stack,
                           flow_field);
    time_t end_time = time(NULL);
  
    vcl_cout << "Level " << L 
             << " took " << difftime(end_time, start_time)
             << " seconds" << vcl_endl;

    // Save intermediate results to disk
    vcl_stringstream filename_ss;
    filename_ss << output_root_ << "flow_map_L" << L << ".png";
    vil_save(flow_field.as_colormap(), filename_ss.str().c_str());

    // Prepare for next level (if there is one).
    if (L > 0)
    {
      flow_field.upscale();
      warp_image_stack(image_stack_pyramid[L-1], flow_field,
                       warped_image_stack);
    }
  }
}
 
//:
void ncm_flow_estimator::set_image_stack(
  vcl_vector< vil_image_view<float> > const& image_stack)
{
  image_stack_ = &image_stack;
  warped_image_stack_ = NULL;
}
 
//:
void ncm_flow_estimator::set_image_stacks(
  vcl_vector< vil_image_view<float> > const& image_stack,
  vcl_vector< vil_image_view<float> > const& warped_image_stack)
{
  image_stack_ = &image_stack;
  warped_image_stack_ = &warped_image_stack;
}
 
//:
void ncm_flow_estimator::set_robustifier(
  ncm_flow_robustifier const& robustifier)
{
  robustifier_ = &robustifier;
}
 
//
//  Private methods
//
 
//:
void ncm_flow_estimator::create_image_pyramid_from(
  // Inputs
  vcl_vector< vil_image_view<vxl_byte> > const& image_stack,
  unsigned n_levels,
  // Outputs
  vcl_vector< vcl_vector< vil_image_view<float> > >& image_stack_pyramid)
{
  const unsigned n_images = image_stack.size();

  image_stack_pyramid.clear();
  image_stack_pyramid.resize(n_levels);

  // Create image pyramids and copy into pyramid structure.
  for (unsigned i = 0; i < n_images; ++i)
  {
    // Allocate space for float_image inside the loop so that a new block of
    // memory is allocated for every image.
    vil_image_view<float> float_image;
    vil_convert_cast(image_stack[i], float_image);

    // Store the original image.
    image_stack_pyramid[0].push_back(float_image);

    // Subsample if needed.
    if (n_levels > 1)
    {
      vil_image_view_base_sptr image_sptr = 
        new vil_image_view<float>(float_image);

      vil_pyramid_image_view<float> pyr(image_sptr, n_levels);

      for (unsigned L = 1; L < n_levels; ++L)
      {
        // Allocate space for subsampled_image inside the loop so that a new
        // block of memory is allocated for every image.
        vil_image_view<float> subsampled_image;
        subsampled_image.deep_copy(pyr(L));
        image_stack_pyramid[L].push_back(subsampled_image);
      }
    }
  }
}
 
//: Interpolate image stack based on given flow field.
//  Note: frame f of warped_image_stack is equivalent to frame f of image_stack
//  projected *back* in time. Therefore, frame f of image_stack should be 
//  closer to frame f+1 of warped_image_stack after the warping.
//  This means the same process can be applied to compute the temporal
//  derivatives, regardless of whether the sequence has been warped or not.
void ncm_flow_estimator::warp_image_stack(
  // Inputs
  vcl_vector< vil_image_view<float> > const& image_stack,
  ncm_flow_field const& flow_field,
  // Outputs
  vcl_vector< vil_image_view<float> >& warped_image_stack)
{
  const unsigned ni = image_stack[0].ni();
  const unsigned nj = image_stack[0].nj();
  const unsigned n_images = image_stack.size();

  warped_image_stack.resize(n_images);
  for (unsigned im = 0; im < n_images; ++im)
    warped_image_stack[im].set_size(ni, nj);
  
  for (unsigned i = 0; i < ni; ++i)
  {
    for (unsigned j = 0; j < nj; ++j)
    {
      double x = static_cast<double>(i) + flow_field.u(i,j);
      if (x < 0)
        x = 0;
      else if (x > ni-1)
        x = ni-1;

      double y = static_cast<double>(j) + flow_field.v(i,j);
      if (y < 0)
        y = 0;
      else if (y > nj-1)
        y = nj-1;

      for (unsigned im = 0; im < n_images; ++im)
      {
        warped_image_stack[im](i,j) = vil_bilin_interp(image_stack[im], x, y);
      }
    }
  }
}
 
//: 
void ncm_flow_estimator::estimate_flow_at_level(
  vcl_vector< vil_image_view<float> > const& image_stack,
  vcl_vector< vil_image_view<float> > const& warped_image_stack,
  ncm_flow_field& flow_field)
{
  // Create the observation mask
  
  // Initialize the variables
  //flow_field.set_orientation();
  //flow_field.set_presence();
  //flow_field.set_displacements();
  
  // Minimize the error function
  ncm_flow_null_robustifier robustifier;
  set_robustifier(robustifier);

  if (warped_image_stack.size() == 0)
  {
    set_image_stack(image_stack);
    estimate(flow_field, /* update = */ false);
  }
  else
  {
    set_image_stacks(image_stack, warped_image_stack);
    estimate(flow_field, /* update = */ true);
  }
}
