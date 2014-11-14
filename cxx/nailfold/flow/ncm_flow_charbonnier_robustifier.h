#ifndef __NCM_FLOW_CHARBONNIER_ROBUSTIFIER_H__
#define __NCM_FLOW_CHARBONNIER_ROBUSTIFIER_H__

#include <vnl/vnl_vector.h>

#include <nailfold/flow/ncm_flow_robustifier.h>

class ncm_flow_charbonnier_robustifier : 
  public ncm_flow_robustifier
{
public: // methods

  ncm_flow_charbonnier_robustifier()
    : epsilon_(1.0),
      a_(1.0)
  {
  }

  void set_epsilon(double epsilon) 
  { 
    epsilon_ = epsilon;
    epsilon_to_a_ = vcl_pow(epsilon_, a_);
  }

  void set_a(double a) 
  { 
    a_ = a;
    epsilon_to_a_ = vcl_pow(epsilon_, a_);
  }

private: // methods

  virtual
  double functor(double x) const;

private: // variables

  double epsilon_;
  double a_;

  double epsilon_to_a_;
};

inline
double ncm_flow_charbonnier_robustifier::functor(
  double x) const
{
  return pow(x^2 + epsilon_, a_) - epsilon_to_a_;
}

#endif __NCM_FLOW_CHARBONNIER_ROBUSTIFIER_H__
