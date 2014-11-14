#include "vil_dtcwt.h"

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
//  Helper functions
//

double complex_magnitude(vcl_complex<double> z)
{
  return vcl_sqrt(z.real()*z.real() + z.imag()*z.imag());
}

double complex_phase(vcl_complex<double> z)
{
  return vcl_atan2(z.imag(), z.real());
}

//
//  Public member functions
//

//
//: Constructor
vil_dtcwt::vil_dtcwt(vil_image_view<vxl_byte> img,
                     unsigned n_levels /* = 1 */)
: n_levels_(n_levels),
  src_left_(0),
  src_top_(0),
  src_ni_(0),
  src_nj_(0),
  max_filter_size_(0)
{
  // set the filters
  clear_filters();
  set_biort_near_sym_b_bp();
  //set_antonini();
  set_qshift_qshift_b_bp();
  
  //clear_filters(); set_dummy_filters(); // temporary

  set_image(img);

  process();
}

//
//: Specify an image to decompose
void vil_dtcwt::set_image(vil_image_view<vxl_byte> img)
{
  // get size of largest filter
  if (h0o_.size() > max_filter_size_)
    max_filter_size_ = h0o_.size();
  if (h1o_.size() > max_filter_size_)
    max_filter_size_ = h1o_.size();
  if (h0b_.size() > max_filter_size_)
    max_filter_size_ = h0b_.size();

  // allow for 1px border added to images whose dimensions are not multiples
  // of four
  unsigned border_size = max_filter_size_ + n_levels_;

  // create three images in which to store image data
  im0_.set_size(img.ni() + 2*border_size, img.nj() + 2*border_size);
  im0_.fill(0);
  im1_.set_size(img.ni() + 2*border_size, img.nj() + 2*border_size);
  im1_.fill(0);
  work_im_.set_size(img.ni() + 2*border_size, img.nj() + 2*border_size);
  work_im_.fill(0);

  // create a view of im0_, and put a <double> copy of the original image 
  // into it.
  src_left_ = border_size;
  src_top_ = border_size;
  vil_image_view<double> src_view = vil_crop(im0_, src_left_, img.ni(),
                                                   src_top_, img.nj());
  vil_convert_cast(img, src_view);

  // update dimensions of source image (to multiples of 2)
  src_ni_ = ((img.ni()+1) /2)*2;
  src_nj_ = ((img.nj()+1) /2)*2;

  // if any extra rows/cols need adding then copy end rows/cols
  src_view = vil_crop(im0_, src_left_, src_ni_,
                            src_top_, src_nj_);
  if (src_ni_ != img.ni())
    copy_col(src_view, src_ni_-2, src_ni_-1);
  if (src_nj_ != img.nj())
    copy_row(src_view, src_nj_-2, src_nj_-1);
}

const vcl_vector< vil_image_view< vcl_complex<double> > >& 
  vil_dtcwt::tree() const 
{ 
  return tree_; 
}

vcl_complex<double> vil_dtcwt::coefficient_at(unsigned x, unsigned y,
                                              unsigned subband, 
                                              unsigned level) const
{
  vcl_vector<double> weights(2);
  switch (subband)
  {
  case 0: weights[0] = -3.0; weights[1] = -1.0; break;
  case 1: weights[0] = -vcl_sqrt(5.0); weights[1] = -vcl_sqrt(5.0); break;
  case 2: weights[0] = -1.0; weights[1] = -3.0; break;
  case 3: weights[0] = 1.0; weights[1] = -3.0; break;
  case 4: weights[0] = vcl_sqrt(5.0); weights[1] = -vcl_sqrt(5.0); break;
  case 5: weights[0] = 3.0; weights[1] = -1.0; break;
  default: assert(false); break;
  }

  for (unsigned i = 0; i < weights.size(); ++i)
  {
    weights[i] *= 3.141592654 / 2.15;
  }

  vcl_complex<double> z_out = 
      interpolate(vil_plane(tree_[level],subband), y, x, level, weights);

  const vcl_complex<double> imag(0,1);

  // Multiply by -1 if level = 0 and subband in (0,2,3,5)
  double level_multiplier = 1.0;
  if (level == 0)
    level_multiplier *= -1.0;

  switch (subband)
  {
  case 0: z_out *= level_multiplier * imag; break;
  case 1: z_out *= -imag; break;
  case 2: z_out *= level_multiplier * imag; break;
  case 3: z_out *= level_multiplier * -1.0; break;
  case 4: break; // z_out *= 1.0 i.e. do nothing
  case 5: z_out *= level_multiplier * -1.0; break;
  default: assert(false); break;
  }

  return z_out;
}

//
// Protected member functions
//

//
// Private member functions
//

void vil_dtcwt::clear_filters()
{
  g0o_.resize(0); g0a_.resize(0); g0b_.resize(0);
  g1o_.resize(0); g1a_.resize(0); g1b_.resize(0);
  g2o_.resize(0); g2a_.resize(0); g2b_.resize(0);
  h0o_.resize(0); h0a_.resize(0); h0b_.resize(0);
  h1o_.resize(0); h1a_.resize(0); h1b_.resize(0);
  h2o_.resize(0); h2a_.resize(0); h2b_.resize(0);
}

//
//: Given a symmetric filter with the first half of the values filled in,
//  reflect the values to the other side (in place).
//  e.g. [a,b,c,0,0] -> [a,b,c,b,a]
void vil_dtcwt::reflect_filter(vcl_vector<double>& f)
{
  int p1 = 0, p2 = f.size()-1;
  for (; p1 < p2; ++p1, --p2)
    f[p2] = f[p1];
}

void vil_dtcwt::set_biort_near_sym_b_bp()
{
  h0o_.resize(13);
    h0o_[0]  = -0.0017578125000000;
    h0o_[1]  = 0.0000000000000000;
    h0o_[2]  = 0.0222656250000000;
    h0o_[3]  = -0.0468750000000000;
    h0o_[4]  = -0.0482421875000000;
    h0o_[5]  = 0.2968750000000000;
    h0o_[6]  = 0.5554687500000000;
    reflect_filter(h0o_);

  g1o_.resize(13);
    g1o_[0]  = -0.0017578125000000;
    g1o_[1]  = -0.0000000000000000;
    g1o_[2]  = 0.0222656250000000;
    g1o_[3]  = 0.0468750000000000;
    g1o_[4]  = -0.0482421875000000;
    g1o_[5]  = -0.2968750000000000;
    g1o_[6]  = 0.5554687500000000;
    reflect_filter(g1o_);

  h1o_.resize(19);
    h1o_[0]  = -0.0000706263950893;
    h1o_[1]  = 0.0000000000000000;
    h1o_[2]  = 0.0013419015066964;
    h1o_[3]  = -0.0018833705357143;
    h1o_[4]  = -0.0071568080357143;
    h1o_[5]  = 0.0238560267857143;
    h1o_[6]  = 0.0556431361607143;
    h1o_[7]  = -0.0516880580357143;
    h1o_[8]  = -0.2997576032366072;
    h1o_[9]  = 0.5594308035714286;
    reflect_filter(h1o_);

  g0o_.resize(19);
    g0o_[0]  = 0.0000706263950893;
    g0o_[1]  = 0.0000000000000000;
    g0o_[2]  = -0.0013419015066964;
    g0o_[3]  = -0.0018833705357143;
    g0o_[4]  = 0.0071568080357143;
    g0o_[5]  = 0.0238560267857143;
    g0o_[6]  = -0.0556431361607143;
    g0o_[7]  = -0.0516880580357143;
    g0o_[8]  = 0.2997576032366072;
    g0o_[9]  = 0.5594308035714286;
    reflect_filter(g0o_);

  h2o_.resize(19);
    h2o_[0] = -0.0003682500256732;
    h2o_[1] = -0.0006222535855797;
    h2o_[2] = -0.0000781782479826;
    h2o_[3] = 0.0041858208470681;
    h2o_[4] = 0.0081917871788836;
    h2o_[5] = -0.0074232740248026;
    h2o_[6] = -0.0615384268799117;
    h2o_[7] = -0.1481582309116905;
    h2o_[8] = -0.1170763016392158;
    h2o_[9] = 0.6529082158435902;
    reflect_filter(h2o_);

  g2o_.resize(19);
    g2o_[0] = -0.0003682500256732;
    g2o_[1] = -0.0006222535855797; // should this (and the others) be negative?
    g2o_[2] = -0.0000781782479826;
    g2o_[3] = 0.0041858208470681;
    g2o_[4] = 0.0081917871788836;
    g2o_[5] = -0.0074232740248026;
    g2o_[6] = -0.0615384268799117;
    g2o_[7] = -0.1481582309116905;
    g2o_[8] = -0.1170763016392158;
    g2o_[9] = 0.6529082158435902;
    reflect_filter(g2o_);
}

void vil_dtcwt::set_biort_antonini()
{
  h0o_.resize(9);
    h0o_[0] = 0.0267487574108101;
    h0o_[1] = -0.0168641184428747;
    h0o_[2] = -0.0782232665289905;
    h0o_[3] = 0.2668641184428729;
    h0o_[4] = 0.6029490182363593;
    reflect_filter(h0o_);

  h1o_.resize(7);
    h1o_[0] = 0.0456358815571251;
    h1o_[1] = -0.0287717631142493;
    h1o_[2] = -0.2956358815571280;
    h1o_[3] = 0.5575435262285023;
    reflect_filter(h1o_);

  g0o_.resize(7);
    g0o_[0] = -0.0456358815571251;
    g0o_[1] = -0.0287717631142493;
    g0o_[2] = 0.2956358815571280;
    g0o_[3] = 0.5575435262285023;
    reflect_filter(g0o_);

  g1o_.resize(9);
    g1o_[0] = 0.0267487574108101;
    g1o_[1] = 0.0168641184428747;
    g1o_[2] = -0.0782232665289905;
    g1o_[3] = -0.2668641184428729;
    g1o_[4] = 0.6029490182363593;
    reflect_filter(g1o_);
}

void vil_dtcwt::set_qshift_qshift_b_bp()
{
  h0b_.resize(14);
    h0b_[0] = -0.0045568956284755;
    h0b_[1] = -0.0054394759372741;
    h0b_[2] = 0.0170252238815540;
    h0b_[3] = 0.0238253847949203;
    h0b_[4] = -0.1067118046866654;
    h0b_[5] = 0.0118660920337970;
    h0b_[6] = 0.5688104207121227;
    h0b_[7] = 0.7561456438925225;
    h0b_[8] = 0.2752953846688820;
    h0b_[9] = -0.1172038876991153;
    h0b_[10] = -0.0388728012688278;
    h0b_[11] = 0.0346603468448535;
    h0b_[12] = -0.0038832119991585;
    h0b_[13] = 0.0032531427636532;

  g0b_.resize(14);
    g0b_[0] = 0.0032531427636532;
    g0b_[1] = -0.0038832119991585;
    g0b_[2] = 0.0346603468448535;
    g0b_[3] = -0.0388728012688278;
    g0b_[4] = -0.1172038876991153;
    g0b_[5] = 0.2752953846688820;
    g0b_[6] = 0.7561456438925225;
    g0b_[7] = 0.5688104207121227;
    g0b_[8] = 0.0118660920337970;
    g0b_[9] = -0.1067118046866654;
    g0b_[10] = 0.0238253847949203;
    g0b_[11] = 0.0170252238815540;
    g0b_[12] = -0.0054394759372741;
    g0b_[13] = -0.0045568956284755;

  h1b_.resize(14);
    h1b_[0]  = -0.0032531427636532;
    h1b_[1]  = -0.0038832119991585;
    h1b_[2]  = -0.0346603468448535;
    h1b_[3]  = -0.0388728012688278;
    h1b_[4]  = 0.1172038876991153;
    h1b_[5]  = 0.2752953846688820;
    h1b_[6]  = -0.7561456438925225;
    h1b_[7]  = 0.5688104207121227;
    h1b_[8]  = -0.0118660920337970;
    h1b_[9]  = -0.1067118046866654;
    h1b_[10] = -0.0238253847949203;
    h1b_[11] = 0.0170252238815540;
    h1b_[12] = 0.0054394759372741;
    h1b_[13] = -0.0045568956284755;

  g1b_.resize(14);
    g1b_[0]  = -0.0045568956284755;
    g1b_[1]  = 0.0054394759372741;
    g1b_[2]  = 0.0170252238815540;
    g1b_[3]  = -0.0238253847949203;
    g1b_[4]  = -0.1067118046866654;
    g1b_[5]  = -0.0118660920337970;
    g1b_[6]  = 0.5688104207121227;
    g1b_[7]  = -0.7561456438925225;
    g1b_[8]  = 0.2752953846688820;
    g1b_[9]  = 0.1172038876991153;
    g1b_[10] = -0.0388728012688278;
    g1b_[11] = -0.0346603468448535;
    g1b_[12] = -0.0038832119991585;
    g1b_[13] = -0.0032531427636532;

  h2b_.resize(14);
    h2b_[0]  = -0.0027716534934754;
    h2b_[1]  = -0.0004329193033811;
    h2b_[2]  = 0.0210100577283097;
    h2b_[3]  = 0.0614446533755929;
    h2b_[4]  = 0.1732414728674278;
    h2b_[5]  = -0.0447647940175083;
    h2b_[6]  = -0.8381378400904721;
    h2b_[7]  = 0.4367873857803173;
    h2b_[8]  = 0.2626918806166865;
    h2b_[9]  = -0.0076247475815125;
    h2b_[10] = -0.0263685613793659;
    h2b_[11] = -0.0254554351814246;
    h2b_[12] = -0.0095951430541611;
    h2b_[13] = -0.0000243562670333;

  g2b_.resize(14);
    g2b_[0]  = -0.0000243562670333;
    g2b_[1]  = -0.0095951430541611;
    g2b_[2]  = -0.0254554351814246;
    g2b_[3]  = -0.0263685613793659;
    g2b_[4]  = -0.0076247475815125;
    g2b_[5]  = 0.2626918806166865;
    g2b_[6]  = 0.4367873857803173;
    g2b_[7]  = -0.8381378400904721;
    g2b_[8]  = -0.0447647940175083;
    g2b_[9]  = 0.1732414728674278;
    g2b_[10] = 0.0614446533755929;
    g2b_[11] = 0.0210100577283097;
    g2b_[12] = -0.0004329193033811;
    g2b_[13] = -0.0027716534934754;
}

void vil_dtcwt::set_dummy_filters()
{
  double scale = 1.0e-9;
  h0o_.resize(3);
    h0o_[0] =  0.5 * scale;
    h0o_[1] =  0.0 * scale;
    h0o_[2] =  0.5 * scale;
  h1o_.resize(3);
    h1o_[0] = -0.5 * scale;
    h1o_[1] =  0.0 * scale;
    h1o_[2] =  0.5 * scale;
  h2o_.resize(3);
    h2o_[0] =  4.0 * scale;
    h2o_[1] =  0.0 * scale;
    h2o_[2] =  2.0 * scale;

  h0b_.resize(2);
    h0b_[0] =  1.0 * scale;
    h0b_[1] =  1.0 * scale;
  h1b_.resize(2);
    h1b_[0] = -1.0 * scale;
    h1b_[1] =  1.0 * scale;
  h2b_.resize(2);
    h2b_[0] = -2.0 * scale;
    h2b_[1] =  1.0 * scale;
}
//
//: 
void vil_dtcwt::disp_im(vil_image_view<double> img)
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

void vil_dtcwt::disp_im(vil_image_view< vcl_complex<double> > img)
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


void vil_dtcwt::testbench()
{
  src_ = im0_;
  dest_ = im1_;

  //level_ = 0; vcl_vector<double> filt = h0o_;
  level_ = 1; vcl_vector<double> filt = h0b_;
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
  vcl_cout<<vcl_endl;
}

//
//:
void vil_dtcwt::process()
{
  // specify any assumptions made in this function
  assert(n_levels_ > 0);
  assert(src_ni_ > 0);
  assert(src_nj_ > 0);

  // clear any existing data
  tree_.resize(0);

  for (level_ = 0; level_ < n_levels_; ++level_)
  {
    // alternate between im0_ and im1_ as source and destination
    if (level_%2 == 0)
    {
      src_ = im0_;
      dest_ = im1_;
    }
    else
    {
      src_ = im1_;
      dest_ = im0_;
    }

    // view of the filtering output
    vil_image_view<double> dest_view;

    // create a view of the source image (including reflected rows)
    vcl_vector<double> lo_filt, hi_filt, bp_filt;
    if (level_ == 0)
    {
      // set appropriate filters
      lo_filt = h0o_;
      hi_filt = h1o_;
      bp_filt = h2o_;

      // at first level, filtering output is same size as source image
      dest_view = vil_crop(dest_, src_left_, src_ni_,
                                  src_top_, src_nj_);
    }
    else
    {
      // set appropriate filters
      lo_filt = h0b_;
      hi_filt = h1b_;
      bp_filt = h2b_;

      // update bounding box and fill copy first/last rows/cols
      if (src_ni_%4 != 0)
      {
        src_left_--;
        src_ni_ += 2;
      }
      if (src_nj_%4 != 0)
      {
        src_top_--;
        src_nj_ += 2;
      }
      vil_image_view<double> src_view;
      src_view = vil_crop(src_, src_left_, src_ni_, 
                                src_top_, src_nj_);
      copy_row(src_view, 1, 0);
      copy_row(src_view, src_nj_-2, src_nj_-1);

      // at subsequent levels, filtering output is half the size as 
      // source image
      dest_view = vil_crop(dest_, src_left_, src_ni_/2,
                                  src_top_, src_nj_/2);
    }

    // create space for complex image
    vil_image_view< vcl_complex<double> > Yh(dest_view.ni()/2, 
                                             dest_view.nj()/2, 6);

    // reflect rows of src_ here for efficiency
    reflect_rows(src_, src_top_, src_nj_, max_filter_size_);

    if ( !h1o_.empty() )
    {
      filt_cols(src_, work_im_, hi_filt);

      // HiLo
      filt_rows(work_im_, dest_, lo_filt); 
      q2c(dest_view, vil_plane(Yh,0), vil_plane(Yh,5));

      if ( h2o_.empty() )
      {
        // HiHi (highpass pair)
        filt_rows(work_im_, dest_, hi_filt);
        q2c(dest_view, vil_plane(Yh,1), vil_plane(Yh,4));
      }
      else 
      {
        // BaBa (bandpass pair)
        filt_cols(src_, work_im_, bp_filt);
        filt_rows(work_im_, dest_, bp_filt);
        q2c(dest_view, vil_plane(Yh,1), vil_plane(Yh,4));
      }
    }

    filt_cols(src_, work_im_, lo_filt);
    if ( !h1o_.empty() )
    {
      // LoHi
      filt_rows(work_im_, dest_, hi_filt);
      q2c(dest_view, vil_plane(Yh,2), vil_plane(Yh,3));
    }

    // add view to dual tree
    tree_.push_back(Yh);
    // LoLo (lowpass pair)
    // At this point, dest_ is a lowpass filtered version of its former self
    filt_rows(work_im_, dest_, lo_filt);

    if (level_ > 0)
    {
      // account for scaling in decimation
      vil_math_scale_values(dest_, 0.5);

      // update bounding box of source image
      src_ni_ = src_ni_ / 2;
      src_nj_ = src_nj_ / 2;
    }
  }

  // Can't do this - dest_ is a image<double> whereas tree_ is a 
  // vector< image<double> >

  //tree_.push_back(dest_);

}

//
//: Convert from quads in y to complex numbers in z
//  This is the function that actually subsamples by 2 in x and y
void vil_dtcwt::q2c(const vil_image_view<double>& img,
                    vil_image_view< vcl_complex<double> >& p_minus_q,
                    vil_image_view< vcl_complex<double> >& p_plus_q)
{
  // Arrange pixels from the corners of the quads into
  // 2 subimages of alternate real and imag pixels.
  //
  //  a----b
  //  |    |
  //  |    |
  //  c----d
  //
  // p = (a + jb) / sqrt(2)
  // q = (d - jc) / sqrt(2)
  //
  // p+q = ((a+d) + j(b-c)) / sqrt(2)
  // p-q = ((a-d) + j(b+c)) / sqrt(2)

  // Assumptions
  // - destination images have been resized correctly already
  assert( p_minus_q.ni() == img.ni()/2 );
  assert( p_minus_q.nj() == img.nj()/2 );
  assert( p_plus_q.ni() == img.ni()/2 );
  assert( p_plus_q.nj() == img.nj()/2 );

  double root_half = vcl_sqrt(0.5);

  for (unsigned j = 0, src_j = 0; j < p_plus_q.nj(); ++j, src_j+=2)
  {
    for (unsigned i = 0, src_i = 0; i < p_plus_q.ni(); ++i, src_i+=2)
    {
      p_minus_q(i,j).real(root_half * (img(src_i,src_j)-img(src_i+1,src_j+1)));
      p_minus_q(i,j).imag(root_half * (img(src_i+1,src_j)+img(src_i,src_j+1)));

      p_plus_q(i,j).real(root_half * (img(src_i,src_j)+img(src_i+1,src_j+1)));
      p_plus_q(i,j).imag(root_half * (img(src_i+1,src_j)-img(src_i,src_j+1)));
    }
  }
}

//
//: Set a row of image to zero
void vil_dtcwt::copy_row(vil_image_view<double>& image,
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
void vil_dtcwt::copy_col(vil_image_view<double>& image,
                         unsigned src, unsigned dest)
{
  copy_row(vil_transpose(image), src, dest);
}

//
//: In-place reflection of n_reflected rows at top and bottom of image
//  Note that the first and last rows are duplicated
void vil_dtcwt::reflect_rows(vil_image_view<double>& image, 
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
void vil_dtcwt::reflect_columns(vil_image_view<double>& image, 
                                unsigned left_column, unsigned n_src_cols,
                                unsigned n_dest_cols)
{
  reflect_rows(vil_transpose(image), left_column, n_src_cols, n_dest_cols);
}

//
//: Filter columns according to the current level of the tree
void vil_dtcwt::filt_cols(vil_image_view<double>& src,
                          vil_image_view<double>& dest,
                          const vcl_vector<double>& kernel)
{
  // assume image rows have already been reflected
  // - filt_cols is always called on the same image at each level of the tree

  if (level_ == 0)
  {
    // kernel parameters
    int n_before = kernel.size()/2;
    int n_after = (kernel.size()-1)/2;

    // vertical strip encompassing source image
    vil_image_view<double> src_view = 
        vil_crop(src, src_left_, src_ni_, 
                      src_top_-n_before, src_nj_+n_before+n_after);

    // vertical strip encompassing destination image
    vil_image_view<double> dest_view = 
        vil_crop(dest, src_left_, src_ni_, 
                       src_top_-n_before, src_nj_+n_before+n_after);

    // filter columns
    filter_columns(src_view, dest_view, kernel);
  }
  else
  {
    // vertical strip encompassing source image
    vil_image_view<double> src_view = 
        vil_crop(src, src_left_, src_ni_, 
                      src_top_-kernel.size(), src_nj_+2*kernel.size());

    // vertical strip encompassing destination image
    vil_image_view<double> dest_view = 
        vil_crop(dest, src_left_, src_ni_, 
                       src_top_, src_nj_);

    coldfilt(src_view, dest_view, kernel);
  }
}

//
//: Filter rows according to the current level of the tree
void vil_dtcwt::filt_rows(vil_image_view<double>& src,
                          vil_image_view<double>& dest,
                          const vcl_vector<double>& kernel)
{
  // kernel parameters
  if (level_ == 0)
  {
    int n_before = kernel.size()/2;
    int n_after = (kernel.size()-1)/2;

    // reflect columns of work_im_ here
    reflect_columns(src, src_left_, src_ni_, max_filter_size_);

    // horizontal strip encompassing source images
    vil_image_view<double> src_view = 
        vil_crop(src, src_left_-n_before, src_ni_+n_before+n_after, 
                      src_top_, src_nj_);

    // horizontal strip encompassing destination images
    vil_image_view<double> dest_view = 
        vil_crop(dest, src_left_-n_before, src_ni_+n_before+n_after, 
                       src_top_, src_nj_);

    filter_rows(src_view, dest_view, kernel);
  }
  else
  {
    // reflect columns of work_im_ here
    reflect_columns(src, src_left_, src_ni_, kernel.size());

    // horizontal strip encompassing source images
    vil_image_view<double> src_view = 
        vil_crop(src, src_left_-kernel.size(), src_ni_+2*kernel.size(), 
                      src_top_, src_nj_);

    // horizontal strip encompassing destination images
    vil_image_view<double> dest_view = 
        vil_crop(dest, src_left_, src_ni_, 
                       src_top_, src_nj_);

    rowdfilt(src_view, dest_view, kernel);
  }
}


//
//: Filter the rows of image X using filter vector h, without decimation.
//  - assumes that the rows have been reflected already
//  - the result is stored in work_im_
void vil_dtcwt::filter_columns(const vil_image_view<double>& src,
                               vil_image_view<double>& dest, 
                               const vcl_vector<double>& kernel)
{
  filter_rows(vil_transpose(src), vil_transpose(dest), kernel);
}

//
//: Filter the columns of image X using filter vector h, without decimation.
//  - dest is same size as src
//  - This does not apply reflection and is therefore suitable for filtering
//    images that have already been reflected
//  - Since vil_convolve_1d filters rows, we filter the transpose of src and
//    store the result in the transpose of dest
void vil_dtcwt::filter_rows(const vil_image_view<double>& src,
                            vil_image_view<double>& dest, 
                            const vcl_vector<double>& kernel)
{
  // Assumptions
  // - src and dest have been resized correctly already
  assert( src.ni() == dest.ni() );
  assert( src.nj() == dest.nj() ); 

  // kernel parameters
  int n_before = kernel.size()/2;
  int n_after = (kernel.size()-1)/2;

  // filter from source -> destination
  vil_convolve_1d(src, 
                  dest,
                  &kernel[n_before],
                  -n_before, n_after,
                  double(), 
                  vil_convolve_no_extend,
                  vil_convolve_no_extend);
}

//
//: Get n_cols columns offset, offset+n, offset+2n, ... from src
vil_image_view<double> 
vil_dtcwt::every_nth_col_of(const vil_image_view<double>& src,
                            unsigned n, unsigned offset /* = 0 */,
                            int n_cols /* = -1 */)
{
  unsigned n_cols_src = 0;
  if (n_cols == -1)
    n_cols_src = src.ni()-offset;
  else
    n_cols_src = n*(n_cols-1) + 1;

  assert( offset+n_cols_src <= src.ni() );

  return vil_decimate( vil_crop(src,offset,n_cols_src,0,src.nj()) ,n,1);
}

//
//: Get n_rows rows offset, offset+n, offset+2n, ... from src
vil_image_view<double> 
vil_dtcwt::every_nth_row_of(const vil_image_view<double>& src,
                            unsigned n, unsigned offset /* = 0 */,
                            int n_rows /* = -1 */)
{
  unsigned n_rows_src = 0;
  if (n_rows == -1)
    n_rows_src = src.nj()-offset;
  else
    n_rows_src = n*(n_rows-1) + 1;

  assert( offset+n_rows_src <= src.nj() );

  return vil_decimate( vil_crop(src,0,src.ni(),offset,n_rows_src), 1,n);
}

//
//: Filters rows
//  - also downsamples by a factor of two
void vil_dtcwt::rowdfilt(const vil_image_view<double>& src,
                         vil_image_view<double>& dest,
                         const vcl_vector<double>& kernel)
{
  // Assumptions:
  // - kernel must be even in length
  assert((kernel.size()%2) == 0);

  // Build subsampled filter kernels
  unsigned n = kernel.size();
  vcl_vector<double> k_fwd_odd(n/2), k_fwd_even(n/2);
  vcl_vector<double> k_rev_odd(n/2), k_rev_even(n/2);

  double dot_prod = 0.0;
  for (unsigned i = 0, j = 0; i < n; i+=2, ++j)
  {
    k_fwd_odd[j] = kernel[i];
    k_fwd_even[j] = kernel[i+1];
    k_rev_odd[j] = kernel[n-i-1];
    k_rev_even[j] = kernel[n-i-2];

    dot_prod += kernel[i] * kernel[n-i-1];
  }

  unsigned n_cols = 0;

  // get four views of source image: every 4th column with different offsets
  vil_image_view<double> src0, src1, src2, src3;
  n_cols = (src.ni() - 3) / 4;
  src0 = every_nth_col_of(src, 4, 5, n_cols);
  src1 = every_nth_col_of(src, 4, 4, n_cols);
  src2 = every_nth_col_of(src, 4, 3, n_cols);
  src3 = every_nth_col_of(src, 4, 2, n_cols);

  // get views of destination image
  vil_image_view<double> dest_odd, dest_even;
  n_cols = dest.ni() / 2; // number of columns in output image
  dest_odd = every_nth_col_of(dest, 2, 0, n_cols/2);
  dest_even = every_nth_col_of(dest, 2, 1, n_cols/2);

  // Perform filtering on columns of extended matrix

  // get a view of those convolved values that are valid
  vil_image_view<double> work_im2(src0.ni(), src0.nj(), src0.nplanes());
  unsigned i_offset = (k_fwd_odd.size()-1)/2;
  vil_image_view<double> valid = vil_crop(work_im2, i_offset, dest_odd.ni(),
                                                    0, dest_odd.nj());

  if (dot_prod > 0)
  {
    filter_rows(src1, work_im2, k_fwd_odd);
    dest_odd.deep_copy(valid);
    filter_rows(src3, work_im2, k_fwd_even);
    vil_math_image_sum(dest_odd, valid, dest_odd);

    filter_rows(src0, work_im2, k_rev_odd);
    dest_even.deep_copy(valid);
    filter_rows(src2, work_im2, k_rev_even);
    vil_math_image_sum(dest_even, valid, dest_even);
  }
  else
  {
    filter_rows(src1, work_im2, k_fwd_odd);
    dest_even.deep_copy(valid);
    filter_rows(src3, work_im2, k_fwd_even);
    vil_math_image_sum(dest_even, valid, dest_even);

    filter_rows(src0, work_im2, k_rev_odd);
    dest_odd.deep_copy(valid);
    filter_rows(src2, work_im2, k_rev_even);
    vil_math_image_sum(dest_odd, valid, dest_odd);
  }
}

void vil_dtcwt::coldfilt(const vil_image_view<double>& src,
                         vil_image_view<double>& dest,
                         const vcl_vector<double>& kernel)
{
  rowdfilt(vil_transpose(src), vil_transpose(dest), kernel);
}

//
//: Interpolate complex coefficients at a given level of the tree.
vcl_complex<double> vil_dtcwt::interpolate(
    vil_image_view< vcl_complex<double> > z_in,
    unsigned row, unsigned column, unsigned level,
    vcl_vector<double> weights) const
{
  // Precompute complex weights
  vcl_vector< vcl_complex<double> > jw(weights.size());
  for (unsigned i = 0; i < weights.size(); ++i)
  {
    jw[i] = vcl_complex<double>(0, weights[i]);
  }

  // Width of a pixel at the desired level. At level zero, this is 0.5 and is
  // halved for every subsequent level.
  // Remember that values start at 0 in C++ but at 1 in Matlab
  const double pixel_width = vcl_pow(2.0, -static_cast<int>(level+1));

  // Offset of pixel (0,0) with respect to the top/left edge of the image at
  // the requested level of the tree. Since the left-hand edge is always at
  // -0.5 (in level coordinates), you need to add half the width of a pixel 
  // to get its centre.
  const double origin = -0.5 + (0.5 * pixel_width);

  // Coordinates of probe pixel in level coordinates
  double xi = origin + column * pixel_width;
  double yi = origin + row * pixel_width;

  // Compute bounds of the pixels surrounding our probe pixel in this level
  // Round down for lower bound
  // (Should be safe as xi and yi should always be nonintegers)
  int x_min = static_cast<int>(xi);
  int y_min = static_cast<int>(yi);
  // Round up for upper bound
  int x_max = x_min + 1;
  int y_max = y_min + 1;

  // Add a border for interpolation
  //const int border = 0; // bilinear interpolation
  const int border = 1; // bicubic interpolation
  x_min -= border;
  x_max += border;
  y_min -= border;
  y_max += border;

  unsigned n_samples_x = x_max - x_min + 1;
  unsigned n_samples_y = y_max - y_min + 1;

  // Translate xi and yi with respect to top left border
  xi -= x_min;
  yi -= y_min;

  vil_image_view<double> z_real(n_samples_x, n_samples_y);
  vil_image_view<double> z_imag(n_samples_x, n_samples_y);
  vil_image_view<double> z_mag(n_samples_x, n_samples_y);

  // Precompute exponentials
  int x = 0, y = 0;
  vcl_vector< vcl_complex<double> > exp_xjw(n_samples_x);
  for (unsigned i = 0, x = x_min; i < n_samples_x; ++i, ++x)
    exp_xjw[i] = vcl_exp(-(static_cast<double>(x+1)) * jw[1]);

  vcl_vector< vcl_complex<double> > exp_yjw(n_samples_y);
  for (unsigned j = 0, y = y_min; j < n_samples_y; ++j, ++y)
    exp_yjw[j] = vcl_exp(-(static_cast<double>(y+1)) * jw[0]);

  for (unsigned i = 0, x = x_min; i < n_samples_x; ++i, ++x)
  {
    int y = y_min;
    for (unsigned j = 0, y = y_min; j < n_samples_y; ++j, ++y)
    {
      vcl_complex<double> z_unwrap = z_in(x,y) * exp_xjw[i] * exp_yjw[j];
      z_mag(i,j) = complex_magnitude(z_unwrap);
      z_real(i,j) = z_unwrap.real() / (z_mag(i,j) + 1e-6);
      z_imag(i,j) = z_unwrap.imag() / (z_mag(i,j) + 1e-6);
    }
  }

  //disp_im(z_mag);

  // Interpolate magnitude, real and imaginary parts, and phase
  // Note that interpolated magnitude is not equal to the magnitude of 
  // z_real_int + i*z_imag_int
  double z_mag_int = vil_bilin_interp_safe(z_mag, xi, yi);
  double z_real_int = vil_bilin_interp_safe(z_real, xi, yi);
  double z_imag_int = vil_bilin_interp_safe(z_imag, xi, yi);
  double z_phase_int = 
      complex_phase(vcl_complex<double>(z_real_int, z_imag_int));

  // Convert back to complex values
  vcl_complex<double> z_out = 
      z_mag_int * vcl_exp((xi+1) * jw[1] + (yi+1) * jw[0] + 
                          vcl_complex<double>(0,z_phase_int));

  return z_out;
}

