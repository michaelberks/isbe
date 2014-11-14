#include "ncm_vessel_detection_sampler.h"
#include "ncm_vessel_patch_trainer.h"

ncm_vessel_detection_sampler::ncm_vessel_detection_sampler()
	:	MAXLEN(255),
		num_images_( 200 ),
		samples_per_image_( 1 ),
		num_samples_( 2000 )
{
}

bool ncm_vessel_detection_sampler::process_args(int argc, char* argv[])
{ 
	img_dir_ = "C:/isbe/asymmetry_project/data/retinograms/DRIVE_clean/training/images/png/";
	pts_dir_ = "C:/isbe/asymmetry_project/data/retinograms/DRIVE_clean/training/vxl_rfs/sample_data/";
	data_path_ = "C:/isbe/asymmetry_project/data/retinograms/DRIVE_clean/training/vxl_rfs/sample_data/sample.dat";
	img_base_name_ = "training_";
	pts_base_name_ = "vessel_xy_";

	img_dir_ += img_base_name_;
	pts_dir_ += pts_base_name_;

	num_images_ = 20;
	samples_per_image_ = 1;
	num_samples_ = 2000;
	return true;
}

//
//: Load an image file
bool ncm_vessel_detection_sampler::readInitialPtsFile(const vcl_string& fname, vcl_vector<vgl_point_2d<double>>& vessel_pts)
{
	vcl_string path = fname;
  int num_pts;
	vcl_ifstream * afs = new vcl_ifstream(path.c_str());
  
	if (afs==0 || !(*afs))
	{
		vcl_cerr<<"Couldn't open "<<path<<vcl_endl;
		delete afs;
		return false;
	}
  
	char *label = new char[MAXLEN];
	char *str1 = new char[MAXLEN];
	
	while ((*afs)>>vcl_ws, !(*afs).eof())
	{
		(*afs)>>label;
    
		if ( (label[0]=='/') &&( label[1]=='/') )
		{
			// Comment line, so read to end
			afs->getline(str1,MAXLEN);
			vcl_cout<<"Comment read from file: " << str1 << vcl_endl;
		}
		else if ( strcmp(label,"version:")==0 )
		{
			// Read in value of goat
			afs->getline(str1,MAXLEN);
			vcl_cout<<"Line read from file: " << str1 << vcl_endl;
		}
		else if ( strcmp(label,"n_points:")==0 )
		{
			// Read in value of goat
			(*afs)>>num_pts;
			vessel_pts.resize( num_pts );
			vcl_cout<<"n_pts in file: " << num_pts << vcl_endl;
		}
		else if (strcmp(label,"{")==0)
		{
				double vx;
				double vy;

				for (int i = 0; i < num_pts; i++) {
					(*afs)>>vx;
					(*afs)>>vy;
					vessel_pts[i ] = vgl_point_2d<double> (vx, vy);
					vcl_cout<<"Read pt " << i << " = " << vessel_pts[i ] << vcl_endl;
				}
		}
		else if (strcmp(label,"}")==0)
			break;
	}

	//Tidy up and return
	delete [] label;
	delete [] str1;
  
	delete afs;
	return true;
}

//
// Main function
int ncm_vessel_detection_sampler::main_fun( int argc, char* argv[] )
{
	process_args(argc, argv);
	
	//Create IO object to load images
	/*vsml_uint16 image_2d_io<double> im_io;*/
	vsml_byte_image_2d_io im_io;
	
	//Create datastream to write out samples - write out the total no of samples on the first line
	vsl_b_ofstream bfs( data_path_ );
	vsl_b_write(bfs, 2*num_samples_);

	//Loop through the training images
	for (int i = 1; i <= num_images_; i++) {

		vcl_string img_path = img_dir_;
		vcl_string pts_path = pts_dir_;

		//Set up the complete name for this image and associated points file
		vcl_stringstream im_ss;
		im_ss << vcl_setw( 2 ) << vcl_setfill( '0' ) << i;
			
		img_path += im_ss.str();
		pts_path += im_ss.str();

		img_path += ".png";
		pts_path += ".pts";
			
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

		//Load in the initial apex
		vcl_vector<vgl_point_2d<double> > vessel_pts;
		bool pts_loaded = readInitialPtsFile( pts_path, vessel_pts );
		if (!pts_loaded)
		{
			vcl_cout << pts_path << " failed to load!" << vcl_endl;
			continue;
		}
		vcl_cout << pts_path << " loaded ok!" << vcl_endl;

		//Create a patch trainer object
		ncm_vessel_patch_trainer patch_trainer;
		patch_trainer.setSamplesPerImage( samples_per_image_ );

		//Set the trainer's output stream and image, the generate data given ground truth apex
		patch_trainer.setStream( &bfs );
		vcl_cout << "Patch trainer stream set" << vcl_endl;
		patch_trainer.setImage( nailfold_image );
		vcl_cout << "Patch trainer image set" << vcl_endl;
		patch_trainer.generateTrainingData( vessel_pts );
		vcl_cout << "Patch trainer data sampled" << vcl_endl;
											
	}
	bfs.close();
	return 0;
}



