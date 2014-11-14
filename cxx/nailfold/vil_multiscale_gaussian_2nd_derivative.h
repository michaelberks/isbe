#ifndef vil_multiscale_gaussian_2nd_derivative_h_
#define vil_multiscale_gaussian_2nd_derivative_h_

//:
// \file
// \brief Class to compute and represent a multiscale decompostion of an image using separable gaussian 2nd derivatives
//
// \author Mike Berks + Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>

#include <vil/algo/vil_convolve_1d.h>
#include <vil/vil_image_view.h>
#include <vimt/vimt_image_pyramid.h>
#include <vimt/vimt_gaussian_pyramid_builder_2d.h>

class vil_multiscale_gaussian_2nd_derivative
{
public:

	enum FeatureType { Complex = 0,
                     OddEven = 1,
										 Mag = 2,
										 Phase = 3,
										 Conj = 4,
										 Odd = 5,
										 Even = 6
	                   };

  //: Constructors
  vil_multiscale_gaussian_2nd_derivative(vimt_image_2d_of<vxl_byte> img,
            unsigned n_levels = 1, double init_scale = 1.0);

  void set_image(vimt_image_2d_of<vxl_byte> img);

  //: Read-only access to tree structure
  const vcl_vector< vil_image_view< double> >& g2d_responses() const; 
	const vcl_vector< vil_image_view< double> >& h2d_responses() const; 

  vcl_complex<double> coefficient_at(unsigned x, unsigned y,
                                     double theta, unsigned level) const;

  //void set_biorthogonal(const vcl_string& biort_type);

  //void set_qshift(const vcl_string& qshift_type);

protected:

private:

  //: Default constructor
  //  Force user to supply an image at construction
  vil_multiscale_gaussian_2nd_derivative();

  void testbench();

  //: Process the image
  void process();

	//Compute the raw filter responses
	void compute_g_responses();
	void compute_h_responses();

	void compute_g_responses(vil_image_view<double> & responses, unsigned level);
	void compute_h_responses(vil_image_view<double> & responses, unsigned level);

	//:Compute steered responses from raw filter responses
	double steered_g_response(unsigned xi, unsigned yi, unsigned level, double theta);
	double steered_h_response(unsigned xi, unsigned yi, unsigned level, double theta);

	bool sample_response_vector(vcl_vector<double> & response_vec,  unsigned xi, unsigned yi, 
		vcl_vector<unsigned> levels, unsigned num_angles = 6, double base_theta = 0,
		unsigned win_size = 1, bool use_g = 1, bool use_h = 1, FeatureType feat = Conj);

  void disp_im(vil_image_view< vcl_complex<double> > img);
  void disp_im(vil_image_view<double> img);

  void copy_row(vil_image_view<double>& image,
                unsigned src, unsigned dest);
  void copy_col(vil_image_view<double>& image,
                unsigned src, unsigned dest);

  void reflect_rows(vil_image_view<double>& image, 
                    unsigned top_row, unsigned n_src_rows, 
                    unsigned n_dest_rows);
  void reflect_columns(vil_image_view<double>& image, 
                       unsigned top_row, unsigned n_src_cols, 
                       unsigned n_dest_cols);

  vcl_complex<double> interpolate(vil_image_view< vcl_complex<double> > z_in,
                                  unsigned row, unsigned column, unsigned level,
                                  vcl_vector<double> weights) const;

  //: Clear all filters
  void clear_filters();
  void reflect_filter(vcl_vector<double>& f);
  void set_dummy_filters();


  //: Set biorthogonal filter
  void set_hilbert_filters();
	void set_gaussian_filters();

  //
  //  Private variables
  //

	vimt_image_pyramid image_pyramid_;
	vimt_gaussian_pyramid_builder_2d<vxl_byte> pyramid_builder_;

	//: Decomposition
	vcl_vector< vil_image_view<double> > g2d_responses_;
	vcl_vector< vil_image_view<double> > h2d_responses_;

  //: Number of levels of decomposition
  unsigned n_levels_;
  unsigned level_;

	//: Initial scale (i.e. sd) of Gaussian kernel
	double init_scale_;
	int half_width_;

	//: Boundary option
	vil_convolve_boundary_option bo_;

  //: Gaussian filters
  vcl_vector<double> g_, dg_, ddg_;
  vcl_vector<double> h2a_, h2b_, h2c_, h2d_;

  //: Workspace images
  vil_image_view<double> im0_;
  vil_image_view<double> work_im_;
  vil_image_view<double> im1_;

  //: Views of workspace images
  vil_image_view<double> src_, dest_;

  unsigned src_left_, src_top_, src_ni_, src_nj_;
	const double Pi_;
  


  //: Biorthogonal thingy
  vcl_vector<double> biort_;

  //: QShift thingy
  vcl_vector<double> qshift_;
};

#endif // vil_dtcwt_h_