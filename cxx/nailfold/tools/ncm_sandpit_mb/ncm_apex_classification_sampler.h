#ifndef NCM_APEX_CLASSIFICATION_SAMPLER
#define NCM_APEX_CLASSIFICATION_SAMPLER

#include <vcl_cmath.h>
#include <vcl_cstddef.h>
#include <vcl_iomanip.h>
#include <vcl_iosfwd.h>
#include <vcl_iostream.h>
#include <vcl_sstream.h>
#include <vcl_string.h>
#include <vcl_vector.h>

#include <vgl/vgl_point_2d.h>

#include <vimt/vimt_image_2d_of.h>

#include <vsml/vsml_uint16_image_2d_io.h>
#include <vsml/vsml_image_2d_io.h>

#include <mrfr/mrfr_1d_forest.h>
#include <mrfr/mrfr_add_all_loaders.h>
#include <vbst/vbst_add_all_loaders.h>

class ncm_apex_classification_sampler
{

	public:

	ncm_apex_classification_sampler();

	int num_images_;
	int samples_per_image_;
	const int MAXLEN;

	//Random forest object
	mrfr_1d_forest forest_;
	vcl_string forest_path_;
	int num_trees_;

	vcl_ofstream class_text_ostream_;

	vcl_string img_dir_;
	vcl_string pts_dir_;
	vcl_string data_path_;
	vcl_string img_base_name_;

	double apex_x_;
	double apex_y_;
	double class_label_;

	bool readInitialPtsFile(const vcl_string& fname);
	bool process_args(int argc, char* argv[]);
	int test_data();
	int load_forest();
	int sample_data();
	void predict_patch(const  vnl_vector<double> & x);

	int main_fun( int argc, char* argv[] );
	

}; //End of class
#endif //NCM_APEX_CLASSIFICATION_SAMPLER


