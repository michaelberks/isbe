#ifndef __NCM_FLOW_ESTIMATOR_SPARSELM_H__
#define __NCM_FLOW_ESTIMATOR_SPARSELM_H__

#include <nailfold/flow/ncm_flow_estimator.h>

class ncm_flow_field;

class ncm_flow_estimator_sparselm 
: public ncm_flow_estimator
{
public: // methods

  //: Estimate the flow field from the input image stack.
  virtual void estimate(
    ncm_flow_field& flow_field,
    bool update);

protected: // variables

private: // methods

private: // variables

};

#endif __NCM_FLOW_ESTIMATOR_SPARSELM_H__
