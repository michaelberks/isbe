#ifndef __NCM_FLOW_COST_FUNCTION_NAG_H__
#define __NCM_FLOW_COST_FUNCTION_NAG_H__

#include <nailfold/flow/ncm_flow_cost_function.h>

class ncm_flow_robustifier;
class ncm_flow_field;

class ncm_flow_cost_function_nag
: public ncm_flow_cost_function
{
public: // methods

  //: Constructor
  ncm_flow_cost_function_nag(
    vcl_vector< vil_image_view<float> > const* image_stack,
    vcl_vector< vil_image_view<float> > const* warped_image_stack = NULL,
    ncm_flow_robustifier const* robustifier = NULL);

  //: Destructor
  ~ncm_flow_cost_function_nag();

protected: // methods

protected: // variables

private: // methods

private: // variables

};

#endif __NCM_FLOW_COST_FUNCTION_NAG_H__
