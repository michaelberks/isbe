#ifndef __NCM_FLOW_COST_FUNCTION_VXL_H__
#define __NCM_FLOW_COST_FUNCTION_VXL_H__

#include <vil/vil_image_view.h>

#include <vnl/vnl_sparse_lst_sqr_function.h>

class ncm_flow_robustifier;
class ncm_flow_field;

class ncm_flow_cost_function_vxl : 
  public vnl_sparse_lst_sqr_function
{
public: // methods

  //: Constructor
  ncm_flow_cost_function_vxl(
    unsigned int num_a,
    unsigned int num_params_per_a,
    unsigned int num_b,
    unsigned int num_params_per_b,
    unsigned int num_params_c,
    unsigned int num_residuals_per_e,
    UseGradient g = use_gradient,
    UseWeights w = no_weights,
    vcl_vector< vil_image_view<vxl_byte> > const* image_stack = NULL,
    ncm_flow_robustifier const* robustifier = NULL);

  //: Construct vnl_sparse_lst_sqr_function.
  // Assumes A consists of \p num_a parameters each of size \p num_params_per_a
  // Assumes B consists of \p num_b parameters each of size \p num_params_per_b
  // Assumes C consists of \p num_params_c parameters
  // \p xmask is a mask for residual availability.  residual e_ij exists only if mask[i][j]==true
  // Assumes each available residual has size \p num_residuals_per_e
  // The optional argument should be no_gradient if the gradf function has not
  // been implemented.  Default is use_gradient.
  ncm_flow_cost_function_vxl(
    unsigned int num_a,
    unsigned int num_params_per_a,
    unsigned int num_b,
    unsigned int num_params_per_b,
    unsigned int num_params_c,
    const vcl_vector<vcl_vector<bool> >& xmask,
    unsigned int num_residuals_per_e,
    UseGradient g = use_gradient,
    UseWeights w = no_weights,
    vcl_vector< vil_image_view<vxl_byte> > const* image_stack = NULL,
    ncm_flow_robustifier const* robustifier = NULL);

  //: Destructor
  ~ncm_flow_cost_function_vxl();

  //: Compute all residuals.
  //  Given the parameter vectors a, b, and c, compute the vector of residuals f.
  //  f has been sized appropriately before the call.
  //  The default implementation computes f by calling fij for each valid
  //  pair of i and j.  You do not need to overload this method unless you
  //  want to provide a more efficient implementation for your problem.
  //virtual void f(
  //  vnl_vector<double> const& a,
  //  vnl_vector<double> const& b,
  //  vnl_vector<double> const& c,
  //  vnl_vector<double>& f);

  //: Compute the residuals from the ith component of a, the jth component of b.
  //  Given the parameter vectors ai, bj, and c, compute the vector of residuals fij.
  //  fij has been sized appropriately before the call.
  virtual void fij(int i, int j,
                   vnl_vector<double> const& ai,
                   vnl_vector<double> const& bj,
                   vnl_vector<double> const& c,
                   vnl_vector<double>& fij);

protected: // variables

private: // methods

  void differentiate_image_stack(
    vcl_vector< vil_image_view<vxl_byte> > const& I);

  //: Compute the brightness constancy penalties
  unsigned long brightness_constancy(
    ncm_flow_field const& flow_field_to_test,
    vnl_vector<double>& f,
    long f_offset);

  //: Compute the spatial smoothness penalties
  unsigned long spatial_smoothness(
    ncm_flow_field const& flow_field_to_test,
    vnl_vector<double>& f,
    long f_offset);

  unsigned long n_observations() const;

  unsigned long n_residuals() const;

private: // variables

  //: Image derivatives.
  vcl_vector< vil_image_view<float> > Ix_;
  vcl_vector< vil_image_view<float> > Iy_;
  vcl_vector< vil_image_view<float> > It_;

  double alpha_;

  //: Pointer to robustifier.
  ncm_flow_robustifier const* robustifier_;

};

#endif __NCM_FLOW_COST_FUNCTION_VXL_H__
