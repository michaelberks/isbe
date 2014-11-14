#ifndef __NCM_FLOW_NULL_ROBUSTIFIER_H__
#define __NCM_FLOW_NULL_ROBUSTIFIER_H__

#include <vnl/vnl_vector.h>

#include <nailfold/flow/ncm_flow_robustifier.h>

class ncm_flow_null_robustifier : 
  public ncm_flow_robustifier
{
public: // methods

  vnl_vector<double> f(
    const vnl_vector<double>& input_vector) const;

private: // methods

  virtual
  double functor(
    double x) const;

};

inline
double ncm_flow_null_robustifier::functor(
  double x) const
{
  return x;
}

#endif __NCM_FLOW_NULL_ROBUSTIFIER_H__
