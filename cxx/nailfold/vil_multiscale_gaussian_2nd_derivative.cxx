#include "vil_multiscale_gaussian_2nd_derivative.h"

#include <vcl_cmath.h>
#include <vcl_cassert.h>
#include <vcl_iomanip.h>
#include <vcl_complex.h>

#include <vil/vil_bilin_interp.h>
#include <vil/vil_convert.h>
#include <vil/vil_transpose.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>
#include <vil/vil_decimate.h>
#include <vil/algo/vil_convolve_1d.h>

//
//  Public member functions
//

//
//: Constructor
vil_multiscale_gaussian_2nd_derivative::vil_multiscale_gaussian_2nd_derivative(vimt_image_2d_of<vxl_byte> img,
                     unsigned n_levels /* = 1 */, double init_scale /* = 1*/)
: n_levels_(n_levels),
  src_left_(0),
  src_top_(0),
  src_ni_(0),
  src_nj_(0),
	Pi_(3.14159265359)
{
  // set the filters
  clear_filters();
  set_hilbert_filters();
	set_gaussian_filters();
  
  //clear_filters(); set_dummy_filters(); // temporary

  set_image(img);

  process();
}

//
//: Specify an image to decompose
void vil_multiscale_gaussian_2nd_derivative::set_image(vimt_image_2d_of<vxl_byte> img)
{
  pyramid_builder_.build(image_pyramid_, img);
	src_ni_ = img.image().ni();
	src_nj_ = img.image().nj();
}

const vcl_vector< vil_image_view< double> >& 
  vil_multiscale_gaussian_2nd_derivative::g2d_responses() const 
{ 
  return g2d_responses_; 
}

const vcl_vector< vil_image_view< double> >& 
  vil_multiscale_gaussian_2nd_derivative::h2d_responses() const 
{ 
  return h2d_responses_; 
}

//
//:
void vil_multiscale_gaussian_2nd_derivative::compute_g_responses()
{
	for (unsigned level = 0; level < n_levels_; ++level)
  {
    //Apply filters to each level of image pyramid
    vil_image_view< double > gr;
		compute_g_responses(gr, level);

		vil_image_view<double> hr;
		compute_g_responses(hr, level);

    // add view to responses tree
    g2d_responses_.push_back(gr);
		h2d_responses_.push_back(hr);
  }
}

 void vil_multiscale_gaussian_2nd_derivative::compute_g_responses(vil_image_view<double> & responses, unsigned level)
{
  set_gaussian_filters();

	//Get image at this level of pyramid
	const vimt_image_2d_of<vxl_byte> &src_img = static_cast<const vimt_image_2d_of<vxl_byte> &>(image_pyramid_(level));
  
	const vil_image_view<vxl_byte> src = src_img.image();

  // create workspace image for reuse throughout
  vil_image_view<double> work_im(src.ni(), src.nj(), src.nplanes());
  vil_image_view<double> work_im_t = vil_transpose(work_im);

	//Resize g responses
	responses.set_size(src.ni(), src.nj(), 3);

  // filter source image three times
  vil_image_view<double> Ixx = vil_plane(responses, 0);
	vil_image_view<double> Iyy = vil_plane(responses, 1);
	vil_image_view<double> Ixy = vil_plane(responses, 2);

  {
    vil_convolve_1d(src, work_im, 
										&ddg_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Ixx);
    vil_convolve_1d(work_im_t, dest_t,
                    &g_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
  }

  {
    vil_convolve_1d(src, work_im,
                    &g_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Iyy);
    vil_convolve_1d(work_im_t, dest_t,
                    &ddg_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
  }

  {
    vil_convolve_1d(src,work_im,
                    &dg_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Ixy);
    vil_convolve_1d(work_im_t,dest_t,
                    &dg_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
  }
  // Invert response later
  //vil_math_scale_values(Ixy, -1.0);

}   

void vil_multiscale_gaussian_2nd_derivative::compute_h_responses()
{
	for (unsigned level = 0; level < n_levels_; ++level)
  {
    //Apply filters to each level of image pyramid
		vil_image_view< double > hr;
		compute_h_responses(hr, level);

    // add view to responses tree
		h2d_responses_.push_back(hr);
  }
}

void vil_multiscale_gaussian_2nd_derivative::compute_h_responses(vil_image_view<double> & responses, unsigned level)
{
	set_hilbert_filters();

	//Get image at this level of pyramid
	const vimt_image_2d_of<vxl_byte> &src_img = static_cast<const vimt_image_2d_of<vxl_byte> &>(image_pyramid_(level));
  
	const vil_image_view<vxl_byte> src = src_img.image();

  // create workspace image for reuse throughout
  vil_image_view<double> work_im(src.ni(), src.nj(), src.nplanes());
  vil_image_view<double> work_im_t = vil_transpose(work_im);

	//Resize h responses
	responses.set_size(src.ni(), src.nj(), 4);

  // filter source image three times
  vil_image_view<double> Ia = vil_plane(responses, 0);
	vil_image_view<double> Ib = vil_plane(responses, 1);
	vil_image_view<double> Ic = vil_plane(responses, 2);
	vil_image_view<double> Id = vil_plane(responses, 3);

  {
    vil_convolve_1d(src, work_im, 
										&h2a_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Ia);
    vil_convolve_1d(work_im_t, dest_t,
                    &h2d_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
  }

  {
    vil_convolve_1d(src, work_im,
                    &h2b_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Ib);
    vil_convolve_1d(work_im_t, dest_t,
                    &h2c_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
  }

  {
    vil_convolve_1d(src,work_im,
                    &h2c_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Ic);
    vil_convolve_1d(work_im_t,dest_t,
                    &h2b_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
	}

	{
    vil_convolve_1d(src,work_im,
                    &h2d_[half_width_], -int(half_width_), half_width_,
                    double(), bo_, bo_);
    vil_image_view<double> dest_t = vil_transpose(Id);
    vil_convolve_1d(work_im_t,dest_t,
                    &h2a_[half_width_],-int(half_width_),half_width_,
                    double(), bo_, bo_);
	}
}

double vil_multiscale_gaussian_2nd_derivative::steered_g_response(unsigned x, unsigned y, unsigned level, double theta)
{
	//Check boundary limits of xi, yi
	assert(x >= 0 && x < src_ni_ && y >= 0 && y < src_nj_);

	//Scale xi, yi to appropriate pyramid level
	const double pixel_width = vcl_pow(2.0, -static_cast<int>(level));
	double xi = double(x)*pixel_width;
	double yi = double(y)* pixel_width;

	//Interpolate Ixx, Iyy and Ixy responses
	double ixx = vil_bilin_interp_safe(g2d_responses_[level], xi, yi, 0);
	double iyy = vil_bilin_interp_safe(g2d_responses_[level], xi, yi, 0);
	double ixy = vil_bilin_interp_safe(g2d_responses_[level], xi, yi, 0);

	//Coumpute response steered at angle theta
	double cc = vcl_pow(vcl_cos(theta), 2);
	double ss = vcl_pow(vcl_sin(theta), 2);
	double s2 = vcl_sin(2*theta);

	double steered_response = ixx*cc + iyy*ss + ixy*s2;
	return steered_response;
}

double vil_multiscale_gaussian_2nd_derivative::steered_h_response(unsigned x, unsigned y, unsigned level, double theta)
{
	//Check boundary limits of xi, yi
	assert(x >= 0 && x < src_nj_ && y >= 0 && y < src_nj_);

	//Scale xi, yi to appropriate pyramid level
	const double pixel_width = vcl_pow(2.0, -static_cast<int>(level));
	double xi = double(x)*pixel_width;
	double yi = double(y)* pixel_width;

	//Interpolate Ixx, Iyy and Ixy responses
	double iha = vil_bilin_interp_safe(h2d_responses_[level], xi, yi, 0);
	double ihb = vil_bilin_interp_safe(h2d_responses_[level], xi, yi, 0);
	double ihc = vil_bilin_interp_safe(h2d_responses_[level], xi, yi, 0);
	double ihd = vil_bilin_interp_safe(h2d_responses_[level], xi, yi, 0);

	//Coumpute response steered at angle theta
	double ccc = vcl_pow(vcl_cos(theta),3); // cos(theta)^3;
  double sss = vcl_pow(vcl_sin(theta),3); // sin(theta)^3
  double ccs3 = 3 * vcl_pow(vcl_cos(theta),2) * vcl_sin(theta); // 3*sin(theta)*cos(theta)^2
  double ssc3 = 3 * vcl_cos(theta) * vcl_pow(vcl_sin(theta),2); // 3*cos(theta)*sin(theta)^2
    
	double steered_response = iha*sss + ihb*ssc3 + ihc*ccs3 + ihd*ccc;
	return steered_response;

}

bool vil_multiscale_gaussian_2nd_derivative::sample_response_vector(vcl_vector<double> & response_vec, unsigned xi, unsigned yi,
																							vcl_vector<unsigned> levels, unsigned num_angles /* = 6*/, double base_theta /*= 0*/,
																							unsigned win_size /*= 1*/, bool use_g /*= 1*/, bool use_h /*= 1*/, FeatureType feat/*= Conj*/)
{
	//Check what levels we're sampling from
	unsigned num_levels = levels.size();

	//If levels is empty, default to all pyramid levels
	if ( !num_levels )
	{
		num_levels = n_levels_;
		levels.resize(num_levels);
		
		unsigned i = 0;
		for ( vcl_vector<unsigned>::iterator itr = levels.begin(), end = levels.end(); itr != end; ++itr, ++i )
		 *itr = i; 
	}

	double theta_step = (2*Pi_/double(num_angles));

	//Compute how big the vector size will be
	unsigned feat_dims;
	switch (feat)
  {
	
		case Complex:
			feat_dims = 2;
			break;
       
		case OddEven:
			feat_dims = 2;
			break;
	
		case Mag:
			feat_dims = 1;
			break;
										 
		case Phase:
			feat_dims = 1;
			break;
										
		case Conj:
			feat_dims = 2;
			break;

		case Odd:
			feat_dims = 1;
			break;

		case Even:
			feat_dims = 1;
			break;
	}

	unsigned vec_size = feat_dims * win_size * win_size * num_levels * num_angles;

	//Allocate space for response vector
	response_vec.resize(vec_size);

	vcl_vector<double>::iterator itr = response_vec.begin();
	double gr, hr, mag, phase;

	for (unsigned ix = xi-1; ix <= xi+1; ix++) {
		for (unsigned iy = yi-1; iy <= yi+1; iy++) {
			for (unsigned il = 0; ix < num_levels; il++) {

				double theta = base_theta;
				for (unsigned ia = 0; ia < num_angles; ia++) {
					//Increment angle
					theta += theta_step;
					
					//Switch on feature type, so we only sample responses we need
					switch (feat)
					{
	
						case Complex:
							//Sample Gaussian and Hilbert responses and compute magnitude and phase
							gr = steered_h_response(ix, iy, levels[il], theta);
							hr = steered_h_response(ix, iy, levels[il], theta);

							mag = vcl_sqrt(gr*gr + hr*hr);
							phase = vcl_atan2(hr, gr);

							*itr = mag; ++itr;
							*itr = phase; ++itr;
							break;
       
						case OddEven:
							//Sample Gaussian and Hilbert responses
							gr = steered_h_response(ix, iy, levels[il], theta);
							hr = steered_h_response(ix, iy, levels[il], theta);

							*itr = gr; ++itr;
							*itr = hr; ++itr;
							break;
	
						case Mag:
							//Sample Gaussian and Hilbert responses and compute magnitude
							gr = steered_h_response(ix, iy, levels[il], theta);
							hr = steered_h_response(ix, iy, levels[il], theta);

							mag = vcl_sqrt(gr*gr + hr*hr);
							*itr = mag; ++itr;
							break;
										 
						case Phase:
							//Sample Gaussian and Hilbert responses and compute phase
							gr = steered_h_response(ix, iy, levels[il], theta);
							hr = steered_h_response(ix, iy, levels[il], theta);

							phase = vcl_atan2(hr, gr);
							*itr = phase; ++itr;
							break;
										
						case Conj:
							//Sample Gaussian and Hilbert responses and convert to conjugate form
							gr = steered_h_response(ix, iy, levels[il], theta);
							hr = steered_h_response(ix, iy, levels[il], theta);

							if (hr < 0)
								hr *= -1;

							mag = vcl_sqrt(gr*gr + hr*hr);
							phase = vcl_atan2(hr, gr);

							*itr = mag; ++itr;
							*itr = phase; ++itr;
							break;

						case Even:
							//Sample Gaussian responses
							gr = steered_h_response(ix, iy, levels[il], theta);
							*itr = gr; ++itr;
							break;

						case Odd:
							//Sample Hilbert responses
							hr = steered_h_response(ix, iy, levels[il], theta);
							*itr = hr; ++itr;
							break;
					}
				}
			}
		}
	}
  return true;
}

vcl_complex<double> vil_multiscale_gaussian_2nd_derivative::coefficient_at(unsigned x, unsigned y,
                                              double theta, 
                                              unsigned level) const
{


  vcl_complex<double> z_out = vcl_complex<double>(0, 0);

  return z_out;
}

//
// Protected member functions
//

//
// Private member functions
//

void vil_multiscale_gaussian_2nd_derivative::clear_filters()
{
  g_.resize(0); dg_.resize(0); ddg_.resize(0);
	h2a_.resize(0); h2b_.resize(0); h2c_.resize(0); h2d_.resize(0);
}

//
//: Given a symmetric filter with the first half of the values filled in,
//  reflect the values to the other side (in place).
//  e.g. [a,b,c,0,0] -> [a,b,c,b,a]
void vil_multiscale_gaussian_2nd_derivative::reflect_filter(vcl_vector<double>& f)
{
  int p1 = 0, p2 = f.size()-1;
  for (; p1 < p2; ++p1, --p2)
    f[p2] = f[p1];
}

void vil_multiscale_gaussian_2nd_derivative::set_gaussian_filters()
{
	const double two_pi = 2*Pi_;

	// if half_width_ is negative then set it to 5*init_scale_
  if (half_width_ < 0)
    half_width_ = static_cast<int>(5.0 * init_scale_);
    
  const int full_width = 2*half_width_ + 1;

	g_.resize(full_width);
	dg_.resize(full_width);
	ddg_.resize(full_width);

	// Compute gaussian kernels from first principles.
  double k = vcl_pow(two_pi, -0.5);
  double sigmasq = init_scale_*init_scale_;
  for (int x = -half_width_; x <= half_width_; ++x)
  {
    unsigned v = x + half_width_;
    g_[v]   = k * vcl_exp(-0.5 * (x*x) / sigmasq);
    dg_[v]  =  -x / sigmasq * g_[v];
    ddg_[v] = (-1 / sigmasq * g_[v]) - (x / sigmasq * dg_[v]);
  }
}

void vil_multiscale_gaussian_2nd_derivative::set_hilbert_filters()
{
  h2a_.resize(14);
    h2a_[0] = -0.0045568956284755;
    h2a_[1] = -0.0054394759372741;
    h2a_[2] = 0.0170252238815540;
    h2a_[3] = 0.0238253847949203;
    h2a_[4] = -0.1067118046866654;
    h2a_[5] = 0.0118660920337970;
    h2a_[6] = 0.5688104207121227;
    h2a_[7] = 0.7561456438925225;
    h2a_[8] = 0.2752953846688820;
    h2a_[9] = -0.1172038876991153;
    h2a_[10] = -0.0388728012688278;
    h2a_[11] = 0.0346603468448535;
    h2a_[12] = -0.0038832119991585;
    h2a_[13] = 0.0032531427636532;

  h2b_.resize(14);
    h2b_[0] = 0.0032531427636532;
    h2b_[1] = -0.0038832119991585;
    h2b_[2] = 0.0346603468448535;
    h2b_[3] = -0.0388728012688278;
    h2b_[4] = -0.1172038876991153;
    h2b_[5] = 0.2752953846688820;
    h2b_[6] = 0.7561456438925225;
    h2b_[7] = 0.5688104207121227;
    h2b_[8] = 0.0118660920337970;
    h2b_[9] = -0.1067118046866654;
    h2b_[10] = 0.0238253847949203;
    h2b_[11] = 0.0170252238815540;
    h2b_[12] = -0.0054394759372741;
    h2b_[13] = -0.0045568956284755;

  h2c_.resize(14);
    h2c_[0]  = -0.0032531427636532;
    h2c_[1]  = -0.0038832119991585;
    h2c_[2]  = -0.0346603468448535;
    h2c_[3]  = -0.0388728012688278;
    h2c_[4]  = 0.1172038876991153;
    h2c_[5]  = 0.2752953846688820;
    h2c_[6]  = -0.7561456438925225;
    h2c_[7]  = 0.5688104207121227;
    h2c_[8]  = -0.0118660920337970;
    h2c_[9]  = -0.1067118046866654;
    h2c_[10] = -0.0238253847949203;
    h2c_[11] = 0.0170252238815540;
    h2c_[12] = 0.0054394759372741;
    h2c_[13] = -0.0045568956284755;

  h2d_.resize(14);
    h2d_[0]  = -0.0045568956284755;
    h2d_[1]  = 0.0054394759372741;
    h2d_[2]  = 0.0170252238815540;
    h2d_[3]  = -0.0238253847949203;
    h2d_[4]  = -0.1067118046866654;
    h2d_[5]  = -0.0118660920337970;
    h2d_[6]  = 0.5688104207121227;
    h2d_[7]  = -0.7561456438925225;
    h2d_[8]  = 0.2752953846688820;
    h2d_[9]  = 0.1172038876991153;
    h2d_[10] = -0.0388728012688278;
    h2d_[11] = -0.0346603468448535;
    h2d_[12] = -0.0038832119991585;
    h2d_[13] = -0.0032531427636532;
}

void vil_multiscale_gaussian_2nd_derivative::set_dummy_filters()
{
  double scale = 1.0e-9;
  h2a_.resize(3);
    h2a_[0] =  0.5 * scale;
    h2a_[1] =  0.0 * scale;
    h2a_[2] =  0.5 * scale;
  h2b_.resize(3);
    h2b_[0] = -0.5 * scale;
    h2b_[1] =  0.0 * scale;
    h2b_[2] =  0.5 * scale;
  h2c_.resize(3);
    h2c_[0] =  4.0 * scale;
    h2c_[1] =  0.0 * scale;
    h2c_[2] =  2.0 * scale;
	h2d_.resize(3);
    h2d_[0] =  2.0 * scale;
    h2d_[1] =  0.0 * scale;
    h2d_[2] =  4.0 * scale;

  g_.resize(2);
    g_[0] =  1.0 * scale;
    g_[1] =  1.0 * scale;
  dg_.resize(2);
    dg_[0] = -1.0 * scale;
    dg_[1] =  1.0 * scale;
  ddg_.resize(2);
    ddg_[0] = -2.0 * scale;
    ddg_[1] =  1.0 * scale;
}
//
//: 
void vil_multiscale_gaussian_2nd_derivative::disp_im(vil_image_view<double> img)
{
  vcl_cout << vcl_setprecision(2);
  for (unsigned j = 0; j < img.nj(); ++j)
  {
    vcl_cout << j << '\t';
    for (unsigned i = 0; i < img.ni(); ++i)
      vcl_cout << vcl_fixed << img(i,j) << '\t';
    vcl_cout << vcl_endl;
  }
  vcl_cout << vcl_endl;
}

void vil_multiscale_gaussian_2nd_derivative::disp_im(vil_image_view< vcl_complex<double> > img)
{
  vcl_cout << vcl_setprecision(2);
  for (unsigned j = 0; j < img.nj(); ++j)
  {
    for (unsigned i = 0; i < img.ni(); ++i)
      vcl_cout << vcl_fixed << img(i,j) << '\t';
    vcl_cout << vcl_endl;
  }
  vcl_cout << vcl_endl;
}


void vil_multiscale_gaussian_2nd_derivative::testbench()
{
  /*src_ = im0_;
  dest_ = im1_;

  //level_ = 0; vcl_vector<double> filt = h0o_;
  level_ = 1; vcl_vector<double> filt = h2a_;
  reflect_rows(src_, src_top_, src_nj_, filt.size());
  disp_im(src_);
  vcl_cout<<vcl_endl;

  filt_cols(src_, work_im_, filt);
  disp_im(src_);
  disp_im(work_im_);
  disp_im(dest_);
  vcl_cout<<vcl_endl;

  filt_rows(work_im_, dest_, filt);
  disp_im(src_);
  disp_im(work_im_);
  disp_im(dest_);
  vcl_cout<<vcl_endl;*/
}

//
//:
void vil_multiscale_gaussian_2nd_derivative::process()
{
  // specify any assumptions made in this function
  assert(n_levels_ > 0);
  assert(src_ni_ > 0);
  assert(src_nj_ > 0);

  // clear any existing data
  g2d_responses_.resize(0);
	h2d_responses_.resize(0);

	//Compute Gaussian and Hilbert transforms
	compute_g_responses();
	compute_h_responses();
  
}


//
//: Set a row of image to zero
void vil_multiscale_gaussian_2nd_derivative::copy_row(vil_image_view<double>& image,
                         unsigned src, unsigned dest)
{
  // assumptions
  assert( src >= 0 );
  assert( src < image.nj() );
  assert( dest >= 0 );
  assert( dest < image.nj() );

  for (unsigned i = 0; i < image.ni(); ++i)
    image(i,dest) = image(i,src);
}

//
//: Set a column of image to zero
void vil_multiscale_gaussian_2nd_derivative::copy_col(vil_image_view<double>& image,
                         unsigned src, unsigned dest)
{
  copy_row(vil_transpose(image), src, dest);
}

//
//: In-place reflection of n_reflected rows at top and bottom of image
//  Note that the first and last rows are duplicated
void vil_multiscale_gaussian_2nd_derivative::reflect_rows(vil_image_view<double>& image, 
                             unsigned top_row, unsigned n_src_rows,
                             unsigned n_dest_rows)
{
  // assumptions
  const unsigned n_above_top = top_row;
  assert( n_above_top >= n_dest_rows );
  const unsigned n_below_bottom = image.nj() - top_row - n_src_rows;
  assert( n_below_bottom >= n_dest_rows );

  int dest_j = 0;
  int src_j = 0;
  int src_diff = 1; // start by increasing src_j

  // Reflect rows at top
  // When the edge of the existing image is reached, start heading back the
  // other way to generate a reflected periodic signal.
  src_j = top_row;
  dest_j = top_row - 1;
  for (unsigned j = 1; j <= n_dest_rows; ++j, --dest_j, src_j += src_diff)
  {
    for (unsigned i = 0; i < image.ni(); ++i)
    {
      image(i,dest_j) = image(i,src_j);
    }

    if ((j % n_src_rows) == 0)
    {
      // We just copied the edge row of the existing image.

      // Step over the edge (if we want to repeat the end row).
      src_j += src_diff;

      // Turn around to head back the opposite way.
      src_diff *= -1;
    }
  }

  // Reflect rows at bottom.
  // Because the top rows have already been reflected repeatedly, a naive
  // reflection here should be sufficient
  const unsigned bottom_row = top_row + n_src_rows - 1;
  src_j = bottom_row;
  dest_j = bottom_row + 1;
  for (unsigned j = 1; j <= n_dest_rows; ++j, ++dest_j, --src_j)
    for (unsigned i = 0; i < image.ni(); ++i)
      image(i,dest_j) = image(i,src_j);
}

//
//: In-place reflection of n_reflected columns at left and right of image
//  Note that the first and last columns are duplicated
void vil_multiscale_gaussian_2nd_derivative::reflect_columns(vil_image_view<double>& image, 
                                unsigned left_column, unsigned n_src_cols,
                                unsigned n_dest_cols)
{
  reflect_rows(vil_transpose(image), left_column, n_src_cols, n_dest_cols);
}



//
//: Interpolate complex coefficients at a given level of the tree.
vcl_complex<double> vil_multiscale_gaussian_2nd_derivative::interpolate(
    vil_image_view< vcl_complex<double> > z_in,
    unsigned row, unsigned column, unsigned level,
    vcl_vector<double> weights) const
{
  
  // Convert back to complex values
  vcl_complex<double> z_out = vcl_complex<double>(1,1);

  return z_out;
}

