#ifndef NCM_BASE_IIH_PATCH_EXTRACTOR_H
#define NCM_BASE_IIH_PATCH_EXTRACTOR_H

//:
// \file
// \Extract vectorized integral image of patch at some given location in an mage
// \author Michael Berks  
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

class ncm_base_iih_patch_extractor
{

    public:

			struct sampling_box
        {
            //: First Point (bottom left of box)
            vgl_point_2d<double> bottom_left_xy;
            
						//: Number of x points sampled (actually sampled along axis_)
            unsigned nx;
            
						//: Number of y points sample (actually sampled normal to axis_)
            unsigned ny;

						//:Direction of x sampling
						vgl_vector_2d<double> du;

						//:Direction of y sampling
						vgl_vector_2d<double> dv;

						//Dimensions of box, relative to box centre
						double yu;
						double yd;
						double xl;
						double xr;

						//Sample spacing
						double spacing;
        };

        
      //construct one
      ncm_base_iih_patch_extractor();

			//Update the dimensions of the sampling box
			void setBoxSizeUp( double d ){ box_.yu = d; };
			void setBoxSizeDown( double d ){ box_.yd = d; };
			void setBoxSizeLeft( double d ){ box_.xl = d; };
			void setBoxSizeRight( double d ){ box_.xr = d; };
			void setBoxSpacing( double s ){ box_.spacing = s; };

			//Set the smoothing parameter
			void setSigma(double smoothSD) {sigma_ = smoothSD;}

			//Update the postion of the sampling box
			void updateBox( const vgl_point_2d<double>& centre, const vgl_vector_2d<double>& axis );

			//Extract patch and copy into vector
			void extractPatch(vnl_vector<double > &patch_vector);

			//Set image
			void setImage(const vimt_image_2d_of<vxl_byte>& img);

			//
			void writeBox(const vcl_string& filename);
      void readBox(const vcl_string& filename);
        
  protected:

		typedef vcl_pair<int,int> ipair_t; //For structuring element offsets as a single object
		vil_structuring_element  convto_vil_structuring_element(  const vcl_vector<ipair_t >& offsets);

		//Image
		vimt_image_2d_of<double > image_;

		//Sampling box
		ncm_base_iih_patch_extractor::sampling_box box_;

		//Smoothing parameter
		double sigma_;
                
};

#endif NCM_BASE_IIH_PATCH_EXTRACTOR_H



