#ifndef __NCM_FLOW_ESTIMATOR_H__
#define __NCM_FLOW_ESTIMATOR_H__

#include <vcl_vector.h>

#include <vil/vil_image_view.h>

class ncm_flow_field;
class ncm_flow_robustifier;

class ncm_flow_estimator
{
public: // methods

  //: Constructor.
  ncm_flow_estimator();

  void set_output_root(
    vcl_string output_root);

  //: Estimate flow from an image stack with n_levels levels.
  void estimate_flow_from(
    // Inputs
    vcl_vector< vil_image_view<vxl_byte> > const& image_stack,
    // Outputs
    ncm_flow_field& flow_field,
    // Optional inputs
    unsigned n_levels /* = 1 */);

  void set_image_stack(
    vcl_vector< vil_image_view<float> > const& image_stack);

  void set_image_stacks(
    vcl_vector< vil_image_view<float> > const& image_stack,
    vcl_vector< vil_image_view<float> > const& warped_image_stack);

  void set_robustifier(
    ncm_flow_robustifier const& robustifier);

  //: Estimate the flow field from the input image stack.
  virtual void estimate(
    ncm_flow_field& flow_field,
    bool update) = 0;

  void estimate_from_pyramid();

protected: // methods

protected: // variables

  //: Pointer to vector of input images.
  vcl_vector< vil_image_view<float> > const* image_stack_;

  //: Pointer to vector of warped input images.
  vcl_vector< vil_image_view<float> > const* warped_image_stack_;

  //: Pointer to robustifier.
  ncm_flow_robustifier const* robustifier_;

private: // methods

  //: 
  void ncm_flow_estimator::create_image_pyramid_from(
    // Inputs
    vcl_vector< vil_image_view<vxl_byte> > const& image_stack,
    unsigned n_levels,
    // Outputs
    vcl_vector< vcl_vector< vil_image_view<float> > >& image_stack_pyramid);

  //: Interpolate image stack based on given flow field.
  //  Note: frame f of warped_image_stack is equivalent to frame f of image_stack
  //  projected *back* in time. Therefore, frame f of image_stack should be 
  //  closer to frame f+1 of warped_image_stack after the warping.
  void ncm_flow_estimator::warp_image_stack(
    // Inputs
    vcl_vector< vil_image_view<float> > const& image_stack,
    ncm_flow_field const& flow_field,
    // Outputs
    vcl_vector< vil_image_view<float> >& warped_image_stack);
 
  //: 
  void estimate_flow_at_level(
    vcl_vector< vil_image_view<float> > const& image_stack,
    vcl_vector< vil_image_view<float> > const& warped_image_stack,
    ncm_flow_field& flow_field);

private: // variables

  vcl_string output_root_;

};

#endif __NCM_FLOW_ESTIMATOR_H__
