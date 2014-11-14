#ifndef __NCM_FLOW_GEMANMCCLURE_ROBUSTIFIER_H__
#define __NCM_FLOW_GEMANMCCLURE_ROBUSTIFIER_H__

#include <vnl/vnl_vector.h>

#include <nailfold/flow/ncm_flow_robustifier.h>

class ncm_flow_gemanmcclure_robustifer : 
  public ncm_flow_robustifier
{
public: // methods

  ncm_flow_gemanmcclure_robustifer()
    : mu_(1.0)
  {
  }

  void set_mu(double mu) 
  { 
    mu_squared_ = mu * mu;
  }

private: // methods

  virtual
  double functor(double x) const;

private: // variables

  double mu_squared_;
}

inline
double ncm_flow_gemanmcclure_robustifier::functor(
  double x) const
{
  return x / (x + mu_squared_);
}

#endif __NCM_FLOW_GEMANMCCLURE_ROBUSTIFIER_H__
