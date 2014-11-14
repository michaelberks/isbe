#ifndef __NCM_FLOW_COST_FUNCTION_H__
#define __NCM_FLOW_COST_FUNCTION_H__

#include <vcl_vector.h>

#include <vil/vil_image_view.h>

class ncm_flow_robustifier;
class ncm_flow_field;

class ncm_flow_cost_function
{
public: // methods

  //: Constructor
  ncm_flow_cost_function(
    vcl_vector< vil_image_view<float> > const* image_stack,
    vcl_vector< vil_image_view<float> > const* warped_image_stack = NULL,
    ncm_flow_robustifier const* robustifier = NULL);

  //: Destructor
  ~ncm_flow_cost_function();

  double* parameters();

  double* observations();

  int n_parameters();
  int n_constrained_parameters();
  int n_observations();

  int Jnnz();
  int JtJnnz();

protected: // methods

  void smoothe_x(
    vil_image_view<float> const& in,
    vil_image_view<float>& out);

  void smoothe_y(
    vil_image_view<float> const& in,
    vil_image_view<float>& out);

  void differentiate_spatially(
    vcl_vector< vil_image_view<float> > const& I);

  void differentiate_temporally(
    vcl_vector< vil_image_view<float> > const& I0,
    vcl_vector< vil_image_view<float> > const& I1);

  //: Compute the brightness constancy penalties
  unsigned long brightness_constancy(
    ncm_flow_field const& flow_field_to_test,
    double* f,
    long f_offset);

  //: Compute the spatial smoothness penalties
  unsigned long spatial_smoothness(
    ncm_flow_field const& flow_field_to_test,
    double* f,
    long f_offset);

  unsigned long n_obs_brightness() const;

  unsigned long n_obs_smoothness() const;

protected: // variables

  //: Image derivatives.
  vcl_vector< vil_image_view<float> > Ix_;
  vcl_vector< vil_image_view<float> > Iy_;
  vcl_vector< vil_image_view<float> > It_;

  double alpha_;

  unsigned ni_, nj_, n_images_;

  //: Pointer to robustifier.
  ncm_flow_robustifier const* robustifier_;

private: // methods

private: // variables

};

#endif __NCM_FLOW_COST_FUNCTION_H__
