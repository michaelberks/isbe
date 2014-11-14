#ifndef NCM_VESSEL_DETECTION_SAMPLER
#define NCM_VESSEL_DETECTION_SAMPLER

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

class ncm_vessel_detection_sampler
{

	public:

	ncm_vessel_detection_sampler();

	int num_images_;
	int samples_per_image_;
	unsigned num_samples_;
	const int MAXLEN;

	vcl_string img_dir_;
	vcl_string pts_dir_;
	vcl_string data_path_;
	vcl_string img_base_name_;
	vcl_string pts_base_name_;

	bool readInitialPtsFile(const vcl_string& fname, vcl_vector<vgl_point_2d<double>>& vessel_pts);
	bool process_args(int argc, char* argv[]);
	
	int main_fun( int argc, char* argv[] );
	

}; //End of class
#endif //NCM_VESSEL_DETECTION_SAMPLER


