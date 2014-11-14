#ifndef NCM_TEST_APEX_AAM_H
#define NCM_TEST_APEX_AAM_H

#include <vcl_cmath.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vcl_vector.h>
#include <vnl/vnl_random.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/algo/vil_line_filter.h>
#include <vil/algo/vil_suppress_non_max.h>
#include <vil/algo/vil_suppress_non_plateau.h>
#include <vil/algo/vil_threshold.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_corners.h>
#include <vil/algo/vil_dog_filter_5tap.h>
#include <vil/algo/vil_gauss_reduce.h>

#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <mbl/mbl_index_sort.h>

#include <nailfold/vil_gaussian_derivatives.h>
#include <nailfold/vil_suppress_non_max_dir.h>
#include <nailfold/vil_find_connected_components.h>

#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>

#include <mdpm/mdpm_add_all_loaders.h>

#include <vsml/vsml_add_all_loaders.h>
#include <vsml/vsml_byte_image_2d_io.h>
#include <vsml/vsml_points_io_2d.h>
#include <vsml/vsml_points_2d.h>
#include <vsml/vsml_shape_model_data.h>
#include <vsml/vsml_parts.h>
#include <vsml/vsml_point_states.h>

//#include <vsml/vsml_sample_stats_1d.h>
//#include <vsml/vsml_aligner.h>

#include <vapm/vapm_add_all_loaders.h>
#include <vapm/vapm_mr_app_model.h>
#include <vapm/vapm_app_model_instance.h>

#include <vaam/vaam_add_all_loaders.h>
#include <vaxm/vaxm_mr_active_model.h>
#include <vaxm/vaxm_mr_active_model_instance.h>
#include <vaxm/vaxm_active_model_instance.h>

#include <vimt/vimt_add_all_binary_loaders.h>
#include <vimt/vimt_image_pyramid.h>
#include <vimt/vimt_image_pyramid_builder.h>
#include <vpdfl/vpdfl_add_all_binary_loaders.h>

class ncm_test_apex_aam
{

	public:
		ncm_test_apex_aam();
		int main_fun( int argc, char* argv[] );

	private:

	//: raw image
	vil_image_view<vxl_byte> raw_image_;

	//: Gaussian 2nd derivative strength
	vil_image_view<double> g2d_strength_;

	//: Gaussian 2nd derivative orientation
	vil_image_view<double> g2d_orientation_;

	//: Image containing components labelled by ID
  vil_image_view<int> component_label_image_;

  //: Vector of component sizes
  vcl_vector<unsigned> component_sizes_;

  //: Image containing components labelled by size
  vil_image_view<vxl_byte> component_size_image_;

	//: Shape model data
	vsml_shape_model_data sm_data;

	//: Active appearance model object
	vaxm_mr_active_model mrGlobalActiveModel_;

	//: Instance of active appearance model
	vaxm_mr_active_model_instance mrGlobalActiveInstance_;

	//: Appearance model object
	vapm_mr_app_model mrAppearanceModel_;
	
	//: Parts model object
	vsml_parts parts_;

	//:Image pyramid
	vimt_image_pyramid *imagePyramid_;

	//: Image input/output
	vsml_byte_image_2d_io im_io;

	//: Points input/output
	vsml_points_io_2d opts_io;

	//: Image to search
	vimt_image_2d_of<vxl_byte> local_image_;

	//: 
	vcl_vector< vgl_point_2d<double> > gblPoints_;

	vcl_vector< vgl_point_2d<double> > initialPoints_;

	vcl_vector< vcl_string > initialPtsFiles_;
	vcl_vector< vcl_string > imagePatchFiles_;

	vcl_string initialPtsDir_;
	vcl_string imagePatchDir_;
	vcl_string outputPtsDir_;
	unsigned int num_candidates_;

	int n_points;
	const int MAXLEN;

	

	void print_usage();
	void configure();
	bool loadImageFrom(vcl_string filename);
	bool loadModelFrom(vcl_string smd_path);
	void reset();
	bool readCandidatesFile(const vcl_string& fname);
	bool readInitialPtsFile(const vcl_string& fname);
	
	

}; //End of class
#endif //NCM_TEST_APEX_AAM_H


