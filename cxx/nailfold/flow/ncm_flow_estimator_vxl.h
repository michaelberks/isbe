#ifndef __NCM_FLOW_ESTIMATOR_VXL_H__
#define __NCM_FLOW_ESTIMATOR_VXL_H__

#include <nailfold/flow/ncm_flow_cost_function_vxl.h>

class ncm_flow_field;

class ncm_flow_estimator_vxl
{
public: // methods

  //: Constructor.
  ncm_flow_estimator_vxl();

  void set_image_stack(
    vcl_vector< vil_image_view<float> > const& image_stack);

  void set_image_stacks(
    vcl_vector< vil_image_view<float> > const& image_stack,
    vcl_vector< vil_image_view<float> > const& warped_image_stack);

  void set_robustifier(
    ncm_flow_robustifier const& robustifier);

  //: Estimate the flow field from the input image stack.
  void estimate(
    ncm_flow_field& flow_field,
    bool update);

protected: // variables

private: // methods

private: // variables

  //: Pointer to vector of input images.
  vcl_vector< vil_image_view<float> > const* image_stack_;

  //: Pointer to vector of warped input images.
  vcl_vector< vil_image_view<float> > const* warped_image_stack_;

  //: Pointer to robustifier.
  ncm_flow_robustifier const* robustifier_;

};

#endif __NCM_FLOW_ESTIMATOR_VXL_H__
