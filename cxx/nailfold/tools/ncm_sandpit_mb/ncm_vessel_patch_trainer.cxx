//:
// \file
// \Perturb model from current position. Can be used to train classifiers of fit quality
// \author: Michael Berks  
//=======================================================================
#include "ncm_vessel_patch_trainer.h"

//------------------------
//Note all points are in the index referencing of the single model
//See also mmf_base_iih_patch_trainer which deals with mapping between global and local submodel indices

ncm_vessel_patch_trainer::ncm_vessel_patch_trainer ():
    pbfs_(0), 
		bDeleteStream_(false),
		y_(1)
{
	//Create a random seed and give it a warm up
  mz_random_.reseed();
  unsigned i_pt=100+mz_random_.lrand32(0,100);
  for(unsigned i=0;i<i_pt;++i)
  {
      unsigned dummy = mz_random_.lrand32(0,100);
  }
	//Set deafult file string to use for out file if no stream provided
  out_file_ = "default_model";
  out_file_ += "-perturbed_fit.bps";

	//Set limits for sampling offsets
	x_offset_lim_ = 15;
	y_offset_lim_ = 15;
	theta_lim_ = 1;

	//set box axis
  axis_.set( 1.0, 0.0 );  //default x-axis
  
	//Set number of samples to extract from each image
	num_samples_ = 100;

}

ncm_vessel_patch_trainer::~ncm_vessel_patch_trainer()
{
    if(bDeleteStream_)
        delete pbfs_;

		//vcl_cout << "Y prior to destruction " << y_ << vcl_endl;
}

void ncm_vessel_patch_trainer::setBoxfDistances(
    double upDist, double downDist,
    double leftDist, double rightDist)
{
	patch_extractor_.setBoxSizeDown( downDist );
	patch_extractor_.setBoxSizeUp( upDist );
	patch_extractor_.setBoxSizeLeft( leftDist );
	patch_extractor_.setBoxSizeRight( rightDist );
}

void ncm_vessel_patch_trainer::generateTrainingData(
		const vcl_vector <vgl_point_2d<double>> & vessel_pts_xy)
{
	//Check if we have a stream to write to - otherwise create one
  if(!pbfs_)
  {
      pbfs_ = new vsl_b_ofstream( out_file_.c_str(),vcl_ios::app );
      bDeleteStream_=true;
      if (!(*pbfs_))
      {
          vcl_string msg("Failed to open output file ");
          msg+= out_file_;
          vcl_cerr << msg << vcl_endl;
          vcl_cout << msg << vcl_endl;
      }
  }
  
	//Do random pertubations and extract patches
	for(unsigned i_pt=0; i_pt < vessel_pts_xy.size(); ++i_pt)
	{

		//Update sampling box with next vessel point
		patch_extractor_.updateBox( vessel_pts_xy[ i_pt ], axis_ );

		//Extract the patch
		patch_extractor_.extractPatch( patch_vector_ );

		//Write the batch to the outstream
		writePatch(*pbfs_);

    //Write out class label
		y_[0] = 0;
		if( i_pt < vessel_pts_xy.size()/2 )
			y_[0] = 1;

		writeY(*pbfs_);
		//vcl_cout << "Writing out point " << i_pt << " class: " << y_[0] << vcl_endl;
	}
}    
                                          
void ncm_vessel_patch_trainer::perturb(
    vgl_point_2d<double>& centre, double& theta,
    vgl_vector_2d<double>& offset_xy)
{
		//Randomly sample x and y offsets
    double x_off = mz_random_.drand64(-x_offset_lim_, x_offset_lim_);
    double y_off = mz_random_.drand64(-y_offset_lim_, y_offset_lim_);
    offset_xy.set(x_off,  y_off);
    
		//Compute new centre given offsets
    vgl_vector_2d<double> norm(-axis_.y(), axis_.x());
    centre += x_off*axis_ + y_off*norm;
    
    //Randomly sample an orientation pertubations
    theta = mz_random_.drand64(-theta_lim_,theta_lim_);      

}

void ncm_vessel_patch_trainer::setImage(const vimt_image_2d_of<vxl_byte>& image)
{
	patch_extractor_.setImage( image );
}

void ncm_vessel_patch_trainer::setAxis(const vgl_vector_2d<double> & axis)
{
    axis_ = axis;
}

void ncm_vessel_patch_trainer::readBox(const vcl_string& filename)
{
    patch_extractor_.readBox( filename );
}

void ncm_vessel_patch_trainer::writeBox(const vcl_string& filename)
{
    patch_extractor_.writeBox( filename );
}

void ncm_vessel_patch_trainer::writePatch(vsl_b_ofstream& bfs)
{
    vsl_b_write(bfs, patch_vector_);
}

void ncm_vessel_patch_trainer::writeY(vsl_b_ofstream& bfs)
{
	vcl_cout << "Y to be written out: " << y_ << vcl_endl; 
  vsl_b_write(bfs, y_);
}

