#ifndef __NCM_FLOW_LORENTZIAN_ROBUSTIFIER_H__
#define __NCM_FLOW_LORENTZIAN_ROBUSTIFIER_H__

#include <vnl/vnl_vector.h>

#include <nailfold/flow/ncm_flow_robustifier.h>

class ncm_flow_lorentzian_robustifier : 
  public ncm_flow_robustifier
{
public: // methods

  ncm_flow_lorentzian_robustifer()
    : nu_(1.0)
  {
  }

  vnl_vector<double> f(const vnl_vector<double>& input_vector) const;

  void set_nu(double nu) 
  { 
    nu_ = nu; 
  }

private: // methods

  virtual
  double functor(double x) const;

private: // variables

  double nu_;
};

inline
double ncm_flow_lorentzian_robustifier::functor(
  double x) const
{
  return vcl_log(1 + 0.5 * (x / nu_)^2);
}

#endif __NCM_FLOW_LORENTZIAN_ROBUSTIFIER_H__
