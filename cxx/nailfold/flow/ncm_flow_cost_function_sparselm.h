#ifndef __NCM_FLOW_COST_FUNCTION_SPARSELM_H__
#define __NCM_FLOW_COST_FUNCTION_SPARSELM_H__

#include <nailfold/flow/ncm_flow_cost_function.h>

class ncm_flow_robustifier;
class ncm_flow_field;
struct splm_crsm;

class ncm_flow_cost_function_sparselm
: public ncm_flow_cost_function
{
public: // methods

  //: Constructor
  ncm_flow_cost_function_sparselm(
    vcl_vector< vil_image_view<float> > const* image_stack,
    vcl_vector< vil_image_view<float> > const* warped_image_stack = NULL,
    ncm_flow_robustifier const* robustifier = NULL);

  //: Destructor
  ~ncm_flow_cost_function_sparselm();

  //: Compute all residuals.
  static
  void f(
    double* parameters, 
    double* residuals, 
    int nvars, 
    int nobs, 
    void* adata);

  //: Specify the Jacobian nonzeros
  static
  void crs_jacobian_sparsity_ST(
    double* parameters, 
    struct splm_crsm* jacobian_map, 
    int nvars, 
    int nobs, 
    void* adata);

  static
  void crs_jacobian(
    double* parameters, 
    struct splm_crsm* jacobian_map, 
    int nvars, 
    int nobs, 
    void* adata);

  //static
  //void ccs_jacobian(
  //  double* parameters, 
  //  struct splm_ccsm* jacobian_map, 
  //  int nvars, 
  //  int nobs, 
  //  void* adata);

protected: // methods

protected: // variables

private: // methods

private: // variables

};

#endif __NCM_FLOW_COST_FUNCTION_SPARSELM_H__
