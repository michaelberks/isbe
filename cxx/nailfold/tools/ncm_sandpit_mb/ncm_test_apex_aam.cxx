#include "ncm_test_apex_aam.h"

ncm_test_apex_aam::ncm_test_apex_aam()
	:	MAXLEN(255),
		initialPoints_( 31 )
{
}

void ncm_test_apex_aam::print_usage()
{
  vcl_cout << "ncm_test_apex_aam -m model_file -c candidates_file -o output_file" << vcl_endl;
  vcl_cout << "Loads in model from model file.\n";
  vcl_cout << "Loads in each candidate apex (image patch and initial points) from candidate file.\n";
  vcl_cout << "Fits AAM to each candidate\n";
  vcl_cout << "Saves predicted points to output file."<<vcl_endl;
  vul_arg_display_usage_and_exit();
}

void ncm_test_apex_aam::configure()
{  
	vpdfl_add_all_binary_loaders(); 
	vsml_add_all_loaders(); 
	mdpm_add_all_loaders(); 
	vapm_add_all_loaders();
	vaam_add_all_loaders();   
	vimt_add_all_binary_loaders();
}

//
//: Load an image file
bool ncm_test_apex_aam::loadImageFrom(vcl_string filename)
{
	// load the image
	vil_image_view<vxl_byte> loaded_image;
	loaded_image = vil_load(filename.c_str());

	if (!loaded_image)
		return false;

	// if image is RGB, convert to greyscale
	vil_image_view<vxl_byte> greyscale_image;
	if (loaded_image.nplanes() > 1)
		vil_convert_planes_to_grey(loaded_image,greyscale_image);
	else
		greyscale_image.deep_copy(loaded_image);

	// remove top row of pixels
	//   Nailfold images typically come with a row of black pixels along the top
	raw_image_ = vil_crop(greyscale_image, 0, greyscale_image.ni(), 
																				 1, greyscale_image.nj()-1);

	// compute Gaussian 2nd derivatives
	vil_gaussian_2nd_derivative(raw_image_, 
															g2d_strength_, g2d_orientation_,
															4.0);

	// suppress non-maximal points and label connected components
	vil_image_view<bool> bool_peaks(raw_image_.ni(), raw_image_.nj());
	vil_suppress_non_max_dir(g2d_strength_, g2d_orientation_, 
													 bool_peaks);
	vil_find_connected_components(bool_peaks,
																component_label_image_, component_sizes_);

	return true;
}

//
//: Load an image file
bool ncm_test_apex_aam::loadModelFrom(vcl_string smd_path)
{
	/*------------------------------------------------------------------------------------
	Load in Shape model data from parameter file
	--------------------------------------------------------------------------------------*/
	//Check file exists
	if ((vul_file::exists(smd_path) && !vul_file::is_directory(smd_path)))
	{
			//Attempt to load APM
			if(!sm_data.readFile(smd_path))
			{
				vcl_cout << "Failed to load " << smd_path << vcl_endl;
				return false;
			}

			else
				vcl_cout << "Loaded shape model data from " << smd_path << vcl_endl;
	}
	else
	{
			//Report that file doesn't exist
			vcl_cout << smd_path << " does not exist " << vcl_endl;
			return false;
	}

	/*------------------------------------------------------------------------------------
		Now load the appearance model (filename obtained from the shape model data just read)
	--------------------------------------------------------------------------------------*/
	
	//Set up file path to APM
	vcl_string apm_path = sm_data.modelDir();
	apm_path += sm_data.modelName();
	apm_path += ".apm";

	//Check file exists
	if ((vul_file::exists(apm_path.c_str()) && !vul_file::is_directory(apm_path.c_str())))
	{
			//Attempt to load APM
			if(!mrAppearanceModel_.loadFile(apm_path.c_str())) //i.e. the .apm file
				vcl_cout << "Failed to load " << apm_path << vcl_endl;
			else
				vcl_cout << "Loaded APM from " << apm_path << vcl_endl;
	}
	else
	{
			//Report that file doesn't exist
			vcl_cout << apm_path << "does not exist " << vcl_endl;
	}
    
	/*------------------------------------------------------------------------------------
		Now load the active appearance model (filename obtained from the SMD data just read)
	--------------------------------------------------------------------------------------*/
  
	//Set up path to AAM file
	vcl_string aam_path = sm_data.modelDir();
	aam_path += sm_data.modelName();
	aam_path += ".aam";

	//Check file exists
	if ((vul_file::exists(aam_path.c_str()) && !vul_file::is_directory(aam_path.c_str())))
	{
			//Attempt to load file
			if(!mrGlobalActiveModel_.loadFile(aam_path.c_str())) //i.e. the .aam
				vcl_cout << "Failed to load " << aam_path << vcl_endl;
			else
				vcl_cout << "Loaded AAM from " << aam_path << vcl_endl;
	}
	else
	{
		//Report that file does not exist
		vcl_cout << aam_path << " does not exist " << vcl_endl;
	}

	// link the active model to the required appearance model
	mrGlobalActiveModel_.set_model(mrAppearanceModel_);

	/*------------------------------------------------------------------------------------
			Attempt to load in parts  (probably parts already in binary mrAppearanceModel_)
	------------------------------------------------------------------------------------*/
 
	//Set up path to parts file
	vcl_string parts_path = sm_data.modelDir();
	parts_path += sm_data.partsFile();
	parts_path += ".parts";

	//Check file exists
	if ((vul_file::exists(parts_path.c_str()) && !vul_file::is_directory(parts_path.c_str())))
	{
			//Attempt to load parts
			if(!parts_.loadFile(parts_path.c_str()))
				vcl_cout << "Failed to load " << parts_path << vcl_endl;
			else
				vcl_cout << "Loaded parts from " << parts_path << vcl_endl;
	}
	else
	{
		//Report that file does not exist
		vcl_cout << parts_path << "does not exist " << vcl_endl;
	}

	mrGlobalActiveInstance_.set_model(mrGlobalActiveModel_);
	n_points = mrGlobalActiveInstance_.instance().points().size();
	return true;
}

void ncm_test_apex_aam::reset()
{
	//Set everthing back to default values
}

bool ncm_test_apex_aam::readCandidatesFile(const vcl_string& fname)
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
	
	num_candidates_ = 0;

	vcl_string initialPtsFile;
	vcl_string imagePatchFile;
	vcl_string colon;

	while ((*afs)>>vcl_ws, !(*afs).eof())
	{
		(*afs)>>label;

		if (strcmp(label,"points_dir:")==0) {
			(*afs)>>initialPtsDir_;
			vcl_cout << initialPtsDir_ << vcl_endl;
		}
		else if (strcmp(label,"image_dir:")==0) {
			(*afs)>>imagePatchDir_;
			vcl_cout << imagePatchDir_ << vcl_endl;
		}
		else if (strcmp(label,"{")==0) {
			
			//We've reached the start of the candidates list, so start recording
			//the files
			while (true) {

				(*afs)>>initialPtsFile;
				
				//If '}' we've recahed the end of the list, so break
				if (strcmp(initialPtsFile.c_str(),"}")==0)
					break;
				
				//Otherwise, get the image patch file
				(*afs)>>colon;
				(*afs)>>imagePatchFile;

				//Add both files to the list
				initialPtsFiles_.push_back( initialPtsFile );
				imagePatchFiles_.push_back( imagePatchFile );

				// increment the candidates counts
				num_candidates_++;
			}
			break;
		}
	}
	vcl_cout << "Loaded " << num_candidates_ << " candidate files for testing" << vcl_endl;

	//Tidy up and return
	delete [] label;
	delete [] str1;
  
	delete afs;
	return true;
}

bool ncm_test_apex_aam::readInitialPtsFile(const vcl_string& fname)
{
	reset();
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
				double apex_x;
				double apex_y;

				for (int i = 0; i <= 30; i++) {
					(*afs)>>apex_x;
					(*afs)>>apex_y;
					initialPoints_[i ] = vgl_point_2d<double> (apex_x, apex_y);
				
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
int ncm_test_apex_aam::main_fun( int argc, char* argv[] )
{
	vul_arg<vcl_string> model_path("-m","Model filename");
  vul_arg<vcl_string> candidates_path("-c","Candidates filename");
  vul_arg<vcl_string> output_path("-o","Output filename");
	vul_arg<vcl_string> quality_path("-q","Model qualities filename");
  vul_arg_parse(argc,argv);

	vcl_cout << model_path.value_ << vcl_endl;
	vcl_cout << candidates_path.value_ << vcl_endl;
	vcl_cout << output_path.value_ << vcl_endl;
	vcl_cout << quality_path.value_ << vcl_endl;

	if ( model_path()=="" )
  {
    print_usage();
    return 0;
  }

	outputPtsDir_ = output_path.value_;
	
	vcl_string quality_path_str;

	if ( quality_path() =="" ) 
	{
		quality_path_str = outputPtsDir_;
		quality_path_str += "model_qualities.txt";
	}
	else
		quality_path_str = quality_path.value_;
	
	//Get list of candidates
	readCandidatesFile( candidates_path.value_ );

	//Loads binary loaders
	configure();

	//Attempt to load models
	loadModelFrom( model_path.value_ );

	//Get the pyramid builder from model instance
	const vimt_image_pyramid_builder& pyramidBuilder = mrGlobalActiveModel_.appModel().imagePyrBuilder();
	imagePyramid_ = pyramidBuilder.new_image_pyramid();

	//Create vector of point states
	vcl_vector<int> gblPointStates_(n_points );
	vcl_fill(gblPointStates_.begin(),gblPointStates_.end(),vsml_PtFixed);

	vcl_string full_im_name;
	vcl_string full_opts_name;
	vcl_string full_ipts_name;

	vcl_ofstream model_qualities_file;
	model_qualities_file.open ( quality_path_str.c_str() );

	for (unsigned int i = 0; i < num_candidates_	; i++) {

			full_im_name = imagePatchDir_;
			full_im_name += imagePatchFiles_[ i ];

			full_opts_name = outputPtsDir_;
			full_opts_name += initialPtsFiles_[ i ];

			full_ipts_name = initialPtsDir_;
			full_ipts_name += initialPtsFiles_[ i ];
			
			bool image_loaded = im_io.loadTheImage(local_image_, full_im_name, "");
			if (image_loaded)
			{

				vcl_cout << full_im_name << " loaded ok!" << vcl_endl;
				pyramidBuilder.build(*imagePyramid_, local_image_);

				//Reset the model to make sure we start from the mean shape
				mrGlobalActiveInstance_.instance().reset();

				//Load in the initial apex
				readInitialPtsFile(full_ipts_name);

				//Convert type
				vsml_points_2d initialFixedPoints;
				initialFixedPoints.setPoints2D(initialPoints_);
			  
				//Copy points and states into instance
				mrGlobalActiveInstance_.setPointState(gblPointStates_,initialFixedPoints);

				//Fit shape model to initial points
				mrGlobalActiveInstance_.instance().fitToPoints(initialFixedPoints);

				//Now we've initialised model clear all states for a free image fit
				mrGlobalActiveInstance_.clearPointState();

				//NB Starts from current model parameters
				//If not set will all be zero ie mean shape in model frame and pose questionable (probably pose of 1st image in training set)
				mrGlobalActiveInstance_.start(*imagePyramid_);
				while (mrGlobalActiveInstance_.update(*imagePyramid_))
				{
			      
				}
				vcl_cout << "model successfully converged" << vcl_endl;
				
				//Get points 				
				const vsml_points_2d& modelShapePoints = static_cast<const vsml_points_2d&>(mrGlobalActiveInstance_.instance().points());
				modelShapePoints.getPoints2D(gblPoints_);
	  
				if (opts_io.savePoints(modelShapePoints, full_opts_name))
					vcl_cout<<"Saved points to " << full_opts_name << vcl_endl;
				else
					vcl_cout<<"Failed to save points to " << full_opts_name << vcl_endl;

				
				model_qualities_file << "Apex " << i << ": " << mrGlobalActiveInstance_.instance().quality( local_image_ ) << vcl_endl;
						
			}
			else
			{
				//vcl_cout << full_im_name << " failed to load!" << vcl_endl;
				break;
			}
		//}
	}
	model_qualities_file.close();

	return 0;
}



