#ifndef NCM_BASE_IIH_PATCH_PREDICTOR_H
#define NCM_BASE_IIH_PATCH_PREDICTOR_H

//:
// \file:
// \Predict vessel apex using RF regression forest
// \author: Michael Berks  
//=======================================================================
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

#include <vsl/vsl_binary_loader.h> 

#include <vsml/vsml_uint16_image_2d_io.h>
#include <vsml/vsml_image_2d_io.h>
#include <vil/vil_load.h>

#include <mrfr/mrfr_vec_forest.h>
#include <mrfr/mrfr_add_all_loaders.h>
#include <vbst/vbst_add_all_loaders.h>
#include "ncm_base_iih_patch_extractor.h"

//MB's comments

class ncm_base_iih_patch_predictor
{
  private:

  public:

            
    //constructor/destructor
    ncm_base_iih_patch_predictor();
    ~ncm_base_iih_patch_predictor();

		int num_images_;
		int num_trees_;
		double spacing_;

		unsigned int ni_;
		unsigned int nj_;
		unsigned pt_counter_;
		unsigned votes_file_counter_;

		vcl_string img_dir_;
		vcl_string res_dir_;
		vcl_string forest_path_;
		vcl_string img_base_name_;

		bool readInitialPtsFile(const vcl_string& fname, vcl_vector<vgl_point_2d<double>>& gt_apex);
		void make_votes_path(const vcl_string& votes_path, vcl_string& votes_path_i);
		bool process_args(int argc, char* argv[]);
		bool process_args();
		void update_prediction( double x, double y );
		int predict_saved_data( vcl_string data_path );
		void predict_patch( );

		vil_image_view<double> vote_image_;
		
		//Create a patch extractor object
		ncm_base_iih_patch_extractor patch_extractor_;

		//Vectorised form of integral image of patch
		vnl_vector<double > patch_vector_;

		//Horizontal axis of patch
		vgl_vector_2d<double>  patch_axis_;

		//Random forest object
		mrfr_vec_forest forest_;
		
		int main_fun( int argc, char* argv[] );

		vcl_ofstream predictor_text_ostream_;

  protected:

        
};

#endif //NCM_BASE_IIH_PATCH_PREDICTOR_H




