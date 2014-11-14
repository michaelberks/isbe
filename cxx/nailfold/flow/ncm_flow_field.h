#ifndef __NCM_FLOW_FIELD_H__
#define __NCM_FLOW_FIELD_H__

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vil/vil_image_view.h>

class ncm_flow_field
{
public: // methods

  ncm_flow_field(
    unsigned ni = 0,
    unsigned nj = 0);

  ~ncm_flow_field();

  void set_size(
    unsigned ni, 
    unsigned nj);

  void set_displacements();

  void set_orientation();

  void set_presence();

  //: Return the horizontal flow at all locations.
  vil_image_view<double> const& u() const;

  //: Return the vertical flow at all locations.
  vil_image_view<double> const& v() const;

  //: Return the horizontal flow at (i,j).
  double u(unsigned i, unsigned j) const;

  //: Return the vertical flow at (i,j).
  double v(unsigned i, unsigned j) const;

  //: Return the index in the parameter vector of the horizontal flow at (i,j).
  int index_of_u(unsigned i, unsigned j);

  //: Return the index in the parameter vector of the vertical flow at (i,j).
  int index_of_v(unsigned i, unsigned j);

  vil_image_view<double> magnitude();
  
  vil_image_view<double> direction();

  vil_image_view<vxl_byte> as_colormap();

  unsigned n_parameters();

  //: Set the flow field to the parameters in v.
  void set_from_vector(
    double const* v);

  //: Increment the flow field by the parameters in v.
  void update_from_vector(
    double* v);

  //: Copy the flow field parameters to three vectors a, b and c (used in
  //  nonlinear optimization).
  void get_to_vector(
    double* v);

  void upscale(unsigned factor = 2);

protected: // methods

protected: // variables

private: // methods

  ////: Ensure size is set at creation by making default constructor private.
  //ncm_flow_field();

  //: Return the greatest flow magnitude across the field
  double max_magnitude();

private: // variables

  vil_image_view<double> u_;
  vil_image_view<double> v_;

  //vnl_matrix< vnl_matrix_fixed<2,2,double> > covariance_;
};

#endif __NCM_FLOW_FIELD_H__
