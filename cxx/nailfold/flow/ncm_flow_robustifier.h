#ifndef __NCM_FLOW_ROBUSTIFIER_H__
#define __NCM_FLOW_ROBUSTIFIER_H__

#include <vnl/vnl_vector.h>

class ncm_flow_robustifier
{
public: // methods

  vnl_vector<double> f(const vnl_vector<double>& input_vector) const;

private: // methods

  virtual
  double functor(double x) const = 0;
};

inline
vnl_vector<double> ncm_flow_robustifier::f(
  const vnl_vector<double>& input_vector) const
{
  vnl_vector<double> output_vector(input_vector.size());

  unsigned n_elements = output_vector.size();
  for (unsigned i = 0; i < n_elements; ++i)
    output_vector[i] = functor(input_vector[i]);

  return output_vector;
}

#endif __NCM_FLOW_ROBUSTIFIER_H__
