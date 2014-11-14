#include "ncm_gaussian_patch_extractor.h"

ncm_gaussian_patch_extractor::ncm_gaussian_patch_extractor ()
:	sigma_(1),
  n_levels_(5),
  n_angles_(6),
  theta_base_(0),
  win_sz_(1),
  feat_type_(vil_multiscale_gaussian_2nd_derivative::FeatureType::Conj)
{ 
}

void ncm_gaussian_patch_extractor::extractPatch(vnl_vector<double > &patch_vector)
{

  //Sample patch from image given sampling box
  vimt_image_2d_of<double> patch;
	vimt_resample_bilin(image_, patch, box_.bottom_left_xy, box_.du, box_.dv, box_.nx, box_.ny);

	/*vcl_stringstream im_ss;
	im_ss << box_.bottom_left_xy.x() << '_' << box_.bottom_left_xy.y();
	
	vcl_string img_path = "C:/isbe/nailfold/data/training/patch_samples/";
	img_path += im_ss.str();
	img_path += ".png";

	vil_image_view<vxl_byte> patch_temp;
	vil_convert_stretch_range(patch.image(), patch_temp);
	vil_save( patch_temp, img_path.c_str() );*/
    
  // Compute integral image of patch
  vil_image_view<double > patch_ii;
  vil_math_integral_image(patch.image(), patch_ii);

	/*vcl_stringstream im_ss;
	im_ss << box_.bottom_left_xy.x() << '_' << box_.bottom_left_xy.y();
	
	vcl_string img_path = "C:/isbe/nailfold/data/training/patch_samples/";
	img_path += im_ss.str();
	img_path += ".png";

	vil_image_view<vxl_byte> patch_temp;
	vil_convert_stretch_range(patch_ii, patch_temp);
	vil_save( patch_temp, img_path.c_str() );*/

  //Convert patch integral image into a vector
	unsigned int ni = patch_ii.ni();
	unsigned int nj = patch_ii.nj();

	assert(patch_ii.istep()==1);
	patch_vector.set_size( ni*nj );
	patch_vector.fill(0.0);

	/*
  double* pvCol = patch_vector.data_block();
  double *v = patch_vector.data_block();
  const double* im = patch_ii.top_left_ptr();
	const vcl_ptrdiff_t jstep = patch_ii.jstep();
  const vcl_ptrdiff_t istep = patch_ii.istep();

  assert(istep==1);
  for (unsigned j=0; j < nj; ++j, im += jstep, v+=ni)
    for (unsigned i=0; i < ni; ++i)
      v[i] = im[i];
	*/

	vcl_copy( patch_ii.top_left_ptr(), patch_ii.top_left_ptr() + patch_vector.size(), patch_vector.data_block() );
	//vcl_memcpy( patch_vector.data_block(), patch_ii.top_left_ptr(),  patch_vector.size()*sizeof(double) );


}

void ncm_gaussian_patch_extractor::setImage(const vimt_image_2d_of<vxl_byte>& image)
{
	//Take a copy of the image in VIL format to filter
	vil_image_view<double> img;
	//img.set_size( image_.image().ni(), image_.image().nj() );

  vil_convert_cast(image.image(), img);

	//Median filter the image?
  vil_structuring_element mask;
  {
      //Symmetric square 3x3
      vcl_vector< ipair_t > offsets;
      for(int i=-1; i<=1; ++i)
      {
          for(int j=-1; j<=1; ++j)
          {
              offsets.push_back(ipair_t(i,j));
          }
      }
      mask = (convto_vil_structuring_element(offsets));
  }
	vil_image_view<double > img_medFilt(img.ni(),img.nj());
  img_medFilt.fill(0.0);
	vil_median(img, img_medFilt, mask);
  
  //Gaussian smooth the median filtered image
  vil_image_view<double > img_gaussFilt(img.ni(),img.nj());
  img_gaussFilt.fill(0.0);
	unsigned halfwidth = unsigned( 3.0*sigma_ );
  vil_gauss_filter_2d(img_medFilt, img_gaussFilt, sigma_, halfwidth);


  //Locally normalise the smoothed image
  const double smoothWidthImage = 20.0; //Compute window size
  double decay_k = vcl_exp(-vcl_log(2.0)/smoothWidthImage);

  vil_image_view<double > img_gaussFilt2, local_sum, local_sum2;

	//Compute local sum of image
  vil_exp_filter_2d(img_gaussFilt, local_sum, decay_k, decay_k); 
	
	//Compute local sum of squared image
  vil_math_image_product(img_gaussFilt, img_gaussFilt, img_gaussFilt2);
  vil_exp_filter_2d(img_gaussFilt2, local_sum2, decay_k, decay_k);

	//Preallocate container for locally normalised
	unsigned ni = img.ni();
	unsigned nj = img.nj();
  vil_image_view<double > img_normalised;
  img_normalised.set_size(ni, nj);

	//Subtract local mean and divide by s.d. at each pixel
	//Set minimum cap on the s.d.
	const double min_sd = 5.0;
  const double min_var = min_sd*min_sd;
  for (unsigned i=0; i<ni; ++i)
  {
      for (unsigned j=0; j<nj; ++j)
      {
          double mean = local_sum(i,j);
          double var = local_sum2(i,j) - mean*mean;
          double sd = min_sd;

          if (var > min_var)
						sd = vcl_sqrt(var);
						
					img_normalised(i,j) = (img(i,j)-mean)/sd;
					//img_normalised(i,j) = img(i,j);
      }
  }

	//Set main image to the normalised image
  image_.image() = img_normalised;
	image_.set_world2im( image.world2im() );
	
	/*vil_image_view<vxl_byte> img_temp;
  vil_convert_stretch_range(img_normalised, img_temp);
	vil_save( img_temp, "normalised_nailfold.png" );*/
	

}

vil_structuring_element  ncm_gaussian_patch_extractor::convto_vil_structuring_element(
    const vcl_vector<ipair_t >& offsets)
{
    vcl_vector<int> pi,pj;
    pi.reserve(offsets.size());
    pj.reserve(offsets.size());
    vcl_vector<ipair_t >::const_iterator ptIter=offsets.begin();
    vcl_vector<ipair_t >::const_iterator ptIterEnd=offsets.end();
    while(ptIter != ptIterEnd)
    {
        pi.push_back(ptIter->first);
        pj.push_back(ptIter->second);
        ++ptIter++;
    }
    vil_structuring_element element(pi,pj);

    return element;
    
}

void ncm_gaussian_patch_extractor::readBox(const vcl_string& filename)
{
    vcl_ifstream ifs(filename.c_str());
    if(!ifs)
    {
        vcl_string msg("Cannot load box from file ");
        msg+=filename;
        //throw qfdv_fileSaveException(msg);
    }
    
    sampling_box b;
    double x,y;
    
    ifs>>x>>y;
		b.bottom_left_xy.set(x,y);
    ifs >> b.nx >> b.ny;
    if(!ifs)
    {
        vcl_string msg("Error boxf point data from file ");
        msg+=filename;
        //throw qfdv_fileSaveException(msg);
    }
    box_ = b;
    vcl_cout << "have read in box location in ncm_gaussian_patch_extractor::readBox" << vcl_endl;
    vcl_cout<<"Box coord " << x << "\t" << y << "\t" <<b.nx<< "\t" << b.ny << vcl_endl;
}

void ncm_gaussian_patch_extractor::writeBox(const vcl_string& filename)
{
    vcl_ofstream ofs(filename.c_str());
    if(!ofs)
    {
        vcl_string msg("Cannot save box to file ");
        msg+=filename;
        //throw qfdv_fileSaveException(msg);
    }
    
		ofs<<box_.bottom_left_xy.x()<<"\t"<<box_.bottom_left_xy.y()<<vcl_endl;
    ofs<<box_.nx<<"\t"<<box_.ny<<vcl_endl;
}