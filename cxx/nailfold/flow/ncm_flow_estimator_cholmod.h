#ifndef __NCM_FLOW_ESTIMATOR_CHOLMOD_H__
#define __NCM_FLOW_ESTIMATOR_CHOLMOD_H__

#include <nailfold/flow/ncm_flow_estimator.h>
#include <nailfold/flow/ncm_flow_cost_function_cholmod.h>

#include <cholmod.h>

class ncm_flow_field;
//class cholmod_sparse;
//class cholmod_dense;
//class cholmod_common;

class ncm_flow_estimator_cholmod
: public ncm_flow_estimator
{
public: // methods

  //: Estimate the flow field from the input image stack.
  virtual void estimate(
    ncm_flow_field& flow_field,
    bool update);

protected: // variables

private: // methods

  void minimize(
    ncm_flow_cost_function_cholmod func,
    vnl_vector<double>& parameters);

  void convert_matrix_to_sparse(
    const vnl_sparse_matrix<double> in,
    cholmod_sparse*& out,
    cholmod_common* common);

  void convert_vector_to_dense(
    const vnl_vector<double> in,
    cholmod_dense*& out,
    cholmod_common* common);

  void set_cholmod_options(
    cholmod_common& c);

  bool converged();

  void solve_system(
    // Input
    vnl_sparse_matrix<double> const& JtJ,
    cholmod_dense* Jte_cm,
    // Input/Output
    cholmod_common& c,
    // Output
    vnl_vector<double>& delta_p);


private: // variables

};

#endif __NCM_FLOW_ESTIMATOR_CHOLMOD_H__
