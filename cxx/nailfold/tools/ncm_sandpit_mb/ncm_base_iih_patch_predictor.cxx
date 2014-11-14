//:
// \file
// \
// \author: Michael Berks  
//=======================================================================
#include "ncm_base_iih_patch_predictor.h"

//------------------------
ncm_base_iih_patch_predictor::ncm_base_iih_patch_predictor ()
{

}

ncm_base_iih_patch_predictor::~ncm_base_iih_patch_predictor()
{

}

bool ncm_base_iih_patch_predictor::process_args(int argc, char* argv[])
{ 
	img_dir_ = "C:/isbe/nailfold/data/test/images/";
	res_dir_ = "C:/isbe/nailfold/data/test/predictions/apex_votes/";
	forest_path_ = "C:/isbe/nailfold/models/vessel/apex_location/test.rf";
	img_base_name_ = "nailfold";

	img_dir_ += img_base_name_;
	res_dir_ += img_base_name_;

	patch_axis_.set( 1.0, 0.0 ); //default axis is x-axis

	spacing_ = 2.0;

	num_trees_ = 1;
	num_images_ = 1;
	return true;
}

bool ncm_base_iih_patch_predictor::process_args()
{ 
	img_dir_ = "C:/isbe/nailfold/data/training/images/";
	res_dir_ = "C:/isbe/nailfold/data/training/predictions/apex_votes/";
	forest_path_ = "C:/isbe/nailfold/models/vessel/apex_location/test.rf";
	img_base_name_ = "nailfold";

	img_dir_ += img_base_name_;
	res_dir_ += img_base_name_;

	patch_axis_.set( 1.0, 0.0 ); //default axis is x-axis

	spacing_ = 2.0;

	num_trees_ = 1;
	num_images_ = 6;
	return true;
}

void ncm_base_iih_patch_predictor::update_prediction( double x, double y )
{
	
	//Extract patch at x, y
	vgl_point_2d<double> patch_centre(x, y);

	//Update the sampling box given this centre
	patch_extractor_.updateBox( patch_centre, patch_axis_ );

	//Extract patch
	patch_extractor_.extractPatch( patch_vector_ );
			
	//Copy reference to trees
	const vcl_vector<mrfr_vec_tree>& trees = forest_.tree();

	//Get prediction
	vnl_vector<double> output( 2 );

	//Loop through trees and get predictions
	for (int i = 0; i < num_trees_; i++)
	{		
		//Evaluate patch vector
		trees[ i ].evaluate( patch_vector_, output);

		//Compute row and column pointed to by output
		unsigned ii = unsigned( x - output[ 0 ] );
		unsigned jj = unsigned( y - output[ 1 ] );
		
		//Increment vote image
		if( (ii>=0) && (jj>=0) && (ii<ni_) && (jj<nj_) )
			vote_image_(ii, jj) += 1;

		//Write output to text
		predictor_text_ostream_ << x - output[ 0 ] << " " << y - output[ 1 ] << vcl_endl;				
	}
}

void ncm_base_iih_patch_predictor::predict_patch( )
{

	//Get prediction
	vnl_vector<double> output( 2 );
	vnl_vector<double> y0_all( num_trees_ );
	vnl_vector<double> y1_all( num_trees_ );

	for (int i = 0; i < num_trees_; i++)
	{
		forest_.tree()[ i ].evaluate( patch_vector_, output);
		y0_all[i] = output[0];
		y1_all[i] = output[1];
	}
	vnl_vector<double> y_out( 2 );
	y_out[0] = y0_all.mean();
	y_out[1] = y1_all.mean();
	vcl_cout << "y output = " << y_out << vcl_endl;

}

int ncm_base_iih_patch_predictor::predict_saved_data( vcl_string data_path )
{
	vsl_b_ifstream bfs( data_path );
  if (!bfs)
  {
    vcl_cout << "Unable to open " << data_path << vcl_endl;
    return false;
  }

  unsigned n_data;
  vsl_b_read(bfs,n_data);
  vcl_cout << "Num points in data: " << n_data << vcl_endl;

  vnl_vector<double> y;

  for (unsigned j=0; j<10; ++j)
  {
    vsl_b_read(bfs, patch_vector_);
    vsl_b_read(bfs, y);

		vcl_cout << "---" << vcl_endl;
		vcl_cout << "X in saved data = " << patch_vector_.size() << " " << patch_vector_.mean() << vcl_endl;
		vcl_cout << "y in saved data = " << y << vcl_endl;
		predict_patch();
  }
	return 0;
}

void ncm_base_iih_patch_predictor::make_votes_path(const vcl_string& votes_path, vcl_string& votes_path_i)
{
	votes_path_i = votes_path;
	
	//Set up the complete name for this image and associated points file
	vcl_stringstream im_ss;
	im_ss << vcl_setw( 3 ) << vcl_setfill( '0' ) << votes_file_counter_;

	votes_path_i += "_";
	votes_path_i += im_ss.str();
	votes_path_i += ".txt";
}

//
// Main function
int ncm_base_iih_patch_predictor::main_fun( int argc, char* argv[] )
{
	process_args();


	//Load in the forest
	vbst_add_all_loaders();
	mrfr_add_all_loaders();
	vsl_b_ifstream bfs_in( forest_path_ );
	vsl_b_read( bfs_in, forest_ );
	vcl_cout << "Forest loaded from " << forest_path_ << vcl_endl;
	num_trees_ = forest_.size();

	//Create IO object to load images
	vsml_byte_image_2d_io im_io;

	//Loop through the training images
	for (int i = 1; i <= num_images_; i++) {

		vcl_string img_path = img_dir_;
		vcl_string mask_path = img_dir_;
		vcl_string vote_im_path = res_dir_;
		vcl_string votes_path = res_dir_;

		//Set up the complete name for this image and associated points file
		vcl_stringstream im_ss;
		im_ss << vcl_setw( 3 ) << vcl_setfill( '0' ) << i;
			
		img_path += im_ss.str();
		img_path += ".bmp";

		mask_path += im_ss.str();
		mask_path += "_mask.bmp";

		vote_im_path += im_ss.str();
		vote_im_path += ".png";

		votes_path += im_ss.str();
			
		//Create image object then try and load it
		//vimt_image_2d_of<double> nailfold_image;
		vimt_image_2d_of<vxl_byte> nailfold_image;
		bool img_loaded = im_io.loadTheImage( nailfold_image, img_path, "" );
		if (!img_loaded)
		{
			vcl_cout << img_path << " failed to load!" << vcl_endl;
			continue;
		}
		vcl_cout << img_path << " loaded ok!" << vcl_endl;

		//Create image object
		vimt_image_2d_of<vxl_byte> nailfold_mask;
		bool mask_loaded = im_io.loadTheImage( nailfold_mask, mask_path, "" );
		if (!mask_loaded)
		{
			vcl_cout << mask_path << " failed to load!" << vcl_endl;
			continue;
		}
		vcl_cout << mask_path << " loaded ok!" << vcl_endl;
		

		//Set image to extractor
		patch_extractor_.setImage( nailfold_image );
		vcl_cout << "Patch trainer image set" << vcl_endl;

		//Get size of image
		ni_ = nailfold_image.image().ni();
		nj_ = nailfold_image.image().nj();

		//Set the vote image to zeros
		vote_image_.set_size(ni_, nj_);
		vote_image_.fill( 0 );

		double buff_x = 25;
		double buff_y = 25;

		//Work out how many steps we need to take in the grid
		unsigned num_steps_i = unsigned( (ni_ - 2*buff_x) / spacing_ );
		unsigned num_steps_j = unsigned( (nj_ - 2*buff_y) / spacing_ );

		//Initialise grid points
		double x = buff_x;
		double y;
		
		//predictor_text_ostream_ << "votes =  zeros(" << num_steps_i*num_steps_j*num_trees_ << ",2);" << vcl_endl;

		pt_counter_ = 1;
		votes_file_counter_ = 1;
		vcl_string votes_path_i;
		make_votes_path(votes_path, votes_path_i);

		//Open a stream to save votes for this image
		predictor_text_ostream_.open ( votes_path_i.c_str() );

		//Loop over x and y, sampling a patch and making predictions at each step
		for (unsigned i = 0; i < num_steps_i; i++ )
		{
			y = buff_y;
			for (unsigned j = 0; j < num_steps_j; j++ )
			{
				//Do stuff
				if ( nailfold_mask.image()( unsigned(x), unsigned(y) ) > 0 )
					update_prediction(x , y);

				pt_counter_ += 1;

				if ( pt_counter_ % 10000 == 0 ) {
					//Close the output file stream
					predictor_text_ostream_.close();
					
					//Increment the file counter
					votes_file_counter_ += 1;

					//Make new file stream
					make_votes_path(votes_path, votes_path_i);
					predictor_text_ostream_.open ( votes_path_i.c_str() );
					vcl_cout << "Opened new file " << votes_path_i << vcl_endl;
				}

				//Update y
				y += spacing_;
			}
			x += spacing_;
		}
		vcl_cout << "Patch stuff done" << vcl_endl;

		//Convert the image so we can view visible output
		vil_image_view<vxl_byte> img_temp;
		vil_convert_stretch_range(vote_image_, img_temp);
		
		//Save the output image
		vcl_cout << "Save to " << vote_im_path << vcl_endl;
		vil_save( img_temp, vote_im_path.c_str() );
		
		//Close the output file stream
		predictor_text_ostream_.close();									
	}

	//All done!
	return 0;
}