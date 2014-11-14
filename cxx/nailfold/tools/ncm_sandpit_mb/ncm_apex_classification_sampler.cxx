#include "ncm_apex_classification_sampler.h"
#include "ncm_base_iih_patch_trainer.h"

ncm_apex_classification_sampler::ncm_apex_classification_sampler()
	:	MAXLEN(255),
		num_images_( 100 ),
		samples_per_image_( 100 )
{
}

bool ncm_apex_classification_sampler::process_args(int argc, char* argv[])
{ 
	img_dir_ = "C:/isbe/nailfold/data/apex_detection/embs/training/aligned/candidate_patches/";
	pts_dir_ = "C:/isbe/nailfold/data/apex_detection/embs/training/aligned/candidate_data/";
	data_path_ = "C:/isbe/nailfold/data/apex_detection/embs/training/aligned/class_sample.dat";
	forest_path_ = "C:/isbe/nailfold/models/vessel/apex_classification/test.rf";
	img_base_name_ = "apex";

	img_dir_ += img_base_name_;
	pts_dir_ += img_base_name_;

	num_images_ = 419;
	samples_per_image_ = 10;
	return true;
}

//
//: Load an image file
bool ncm_apex_classification_sampler::readInitialPtsFile(const vcl_string& fname)
{
	vcl_string path = fname;
  
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
		else if (strcmp(label,"version:")==0 || strcmp(label,"n_points:")==0)
		{
			// Read in value of goat
			afs->getline(str1,MAXLEN);
			vcl_cout<<"Line read from file: " << str1 << vcl_endl;
		}
		else if (strcmp(label,"{")==0)
		{
			(*afs)>>apex_x_;
			(*afs)>>apex_y_;
			(*afs)>>class_label_;
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

int ncm_apex_classification_sampler::sample_data( ) {
	//Create IO object to load images
	/*vsml_uint16 image_2d_io<double> im_io;*/
	vsml_byte_image_2d_io im_io;
	
	//Create datastream to write out samples - write out the total no of samples on the first line
	vsl_b_ofstream bfs( data_path_ );
	vsl_b_write(bfs, unsigned( num_images_ * samples_per_image_ ));

	//Loop through the training images
	for (int i = 1; i <= num_images_; i++) {

		vcl_string img_path = img_dir_;
		vcl_string pts_path = pts_dir_;

		//Set up the complete name for this image and associated points file
		vcl_stringstream im_ss;
		im_ss << vcl_setw( 4 ) << vcl_setfill( '0' ) << i;
			
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
		bool pts_loaded = readInitialPtsFile( pts_path );
		if (!pts_loaded)
		{
			vcl_cout << pts_path << " failed to load!" << vcl_endl;
			continue;
		}
		vcl_cout << pts_path << " loaded ok!" << vcl_endl;

		//Create a patch trainer object
		ncm_base_iih_patch_trainer patch_trainer;
		patch_trainer.setSamplesPerImage( samples_per_image_ );

		//Set the trainer's output stream and image, the generate data given ground truth apex
		patch_trainer.setStream( &bfs );
		vcl_cout << "Patch trainer stream set" << vcl_endl;
		patch_trainer.setImage( nailfold_image );
		vcl_cout << "Patch trainer image set" << vcl_endl;
		patch_trainer.generateClassData( vgl_point_2d<double> (apex_x_, apex_y_) , class_label_);
		vcl_cout << "Patch trainer data sampled" << vcl_endl;
								
	}
	bfs.close();
	return 0;
}

int ncm_apex_classification_sampler::test_data( ) {
	
	load_forest();

	vsl_b_ifstream bfs_in( data_path_ );
	if (!bfs_in)
	{
		vcl_cout << "Unable to open " << data_path_ << vcl_endl;
		return false;
	}

	unsigned n_data;
	vsl_b_read(bfs_in,n_data);
	vcl_cout << "Num points in data: " << n_data << vcl_endl;

	vnl_vector<double> x;
	vnl_vector<double> y;

	class_text_ostream_.open("C:/isbe/nailfold/data/apex_detection/embs/test/aligned/class_output.txt");

	for (unsigned j=0; j<n_data; j++)
	{
		vsl_b_read(bfs_in, x);
		vsl_b_read(bfs_in, y);

		/*if (j % samples_per_image_ == 0 ) {
			vcl_cout << "---" << vcl_endl;
			vcl_cout << "X in saved data = " << x.size() << " " << x.mean() << vcl_endl;
			vcl_cout << "y in saved data = " << y << vcl_endl;
		}*/
		vcl_cout << "y in = " << y << " ,";
		
		predict_patch( x );
	}
	class_text_ostream_.close();
	return 0;
}

void ncm_apex_classification_sampler::predict_patch(const  vnl_vector<double> & x)
{

	//Get prediction
	double y_tree;
	vnl_vector<double> y_all( num_trees_ );

	for (int i = 0; i < num_trees_; i++)
	{
		forest_.tree()[ i ].evaluate( x, y_tree);
		y_all[i] = y_tree;
	}
	double y_out = y_all.mean();
	vcl_cout << "y output = " << y_out << vcl_endl;
	class_text_ostream_ << y_out << vcl_endl;

}

int ncm_apex_classification_sampler::load_forest( ) {
	
	//Load in the forest
	vbst_add_all_loaders();
	mrfr_add_all_loaders();
	vsl_b_ifstream bfs_in( forest_path_ );
	vsl_b_read( bfs_in, forest_ );
	vcl_cout << "Forest loaded from " << forest_path_ << vcl_endl;
	num_trees_ = forest_.size();

	return 0;
}


//
// Main function
int ncm_apex_classification_sampler::main_fun( int argc, char* argv[] )
{
	process_args(argc, argv);

	sample_data();
	//test_data();

	
	return 0;
}



