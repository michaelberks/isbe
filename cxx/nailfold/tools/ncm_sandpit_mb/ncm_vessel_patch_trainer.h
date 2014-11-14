#ifndef NCM_VESSEL_PATCH_TRAINER
#define NCM_VESSEL_PATCH_TRAINER

//:
// \file:
// \Perturb shape from current position. Can be used to train classifiers using Brief type features
// \author: Michael Berks  
//=======================================================================
#include <vcl_algorithm.h>
#include <vcl_functional.h>
#include <vcl_numeric.h>
#include <vcl_iterator.h>

#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_2d.h>

#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_plane.h>
#include <vil/vil_math.h>
#include <vil/vil_save.h>
#include <vil/algo/vil_median.h>
#include <vil/algo/vil_exp_filter_2d.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vil/algo/vil_structuring_element.h>

#include <vimt/vimt_image_2d_of.h>
#include <vimt/vimt_resample_bilin.h>

#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector_ref.h>

#include <vsml/vsml_byte_image_2d_io.h>
#include <vsml/vsml_points_io_2d.h>
#include <vsml/vsml_image_2d_io.h> 

#include "ncm_base_iih_patch_extractor.h"

//MB's comments

class ncm_vessel_patch_trainer
{
  private:
    vcl_string out_file_;
    protected:

        //: Number of trainig perturbations of actual shape
        unsigned num_samples_;
        vnl_random mz_random_;

    public:
           
      //constructor/destructor
      ncm_vessel_patch_trainer();
      ~ncm_vessel_patch_trainer();
      
      //
			void generateTrainingData(const vcl_vector < vgl_point_2d<double> > & apex_centre_xy);

      void setImage(const vimt_image_2d_of<vxl_byte>& img);
      
			void setSamplesPerImage(unsigned n) {num_samples_ = n;}
      void setSigma(double smoothSD);
      void setBoxfDistances(double inDist, double outDist, double leftDist, double  rightDist);
      void setAxis(const vgl_vector_2d<double> & axis); //NB axis should be unit vector

      void setStream(vsl_b_ofstream* pbfs) {pbfs_=pbfs;}

      void writePatch(vsl_b_ofstream& bfs);
			void writeY(vsl_b_ofstream& bfs);

			void writeBox(const vcl_string& filename);
      void readBox(const vcl_string& filename);

  protected:

      //: Perturb starting state prior to randomised fit
      void perturb(vgl_point_2d<double>& centre,double& theta,
                   vgl_vector_2d<double>& dp);

			//Object to extract patches from an image
			ncm_base_iih_patch_extractor patch_extractor_;

			//Vectorised form of integral image of patch
			vnl_vector<double > patch_vector_;

			//Random offset from apex centre to be recorded
			vnl_vector<double > y_;
      
      //Unit vector for of unperturbed axis
      vgl_vector_2d<double > axis_;

			//Limits of random pertubtaions
			double x_offset_lim_;
			double y_offset_lim_;
			double theta_lim_;
      
			//Output stream for data
      vsl_b_ofstream* pbfs_;
      bool bDeleteStream_;

};

#endif //NCM_VESSEL_PATCH_TRAINER




