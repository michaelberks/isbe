#ifndef vil_dtcwt_h_
#define vil_dtcwt_h_

//:
// \file
// \brief Class to compute and represent the Dual Tree Complex Wavelet Transform
//        Each level of the tree is stored as a 6-plane image - one plane for
//        each subband. Multiple levels are stored as a vcl_vector of images.
//
//        There could be an argument for using the model-builder design pattern
//        but the current implementation works for now.
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>

#include <vil/vil_image_view.h>

class vil_dtcwt
{
public:

  //: Constructors
  vil_dtcwt(vil_image_view<vxl_byte> img,
            unsigned n_levels = 1);

  void set_image(vil_image_view<vxl_byte> img);

  //: Read-only access to tree structure
  const vcl_vector< vil_image_view< vcl_complex<double> > >& tree() const;

  vcl_complex<double> coefficient_at(unsigned x, unsigned y,
                                     unsigned subband, unsigned level) const;

  //void set_biorthogonal(const vcl_string& biort_type);

  //void set_qshift(const vcl_string& qshift_type);

protected:

private:

  //: Default constructor
  //  Force user to supply an image at construction
  vil_dtcwt();

  void testbench();

  //: Process the image
  void process();

  void disp_im(vil_image_view< vcl_complex<double> > img);
  void disp_im(vil_image_view<double> img);

  void copy_row(vil_image_view<double>& image,
                unsigned src, unsigned dest);
  void copy_col(vil_image_view<double>& image,
                unsigned src, unsigned dest);

  vil_image_view<double> every_nth_row_of(const vil_image_view<double>& src,
                                          unsigned n, unsigned offset = 0,
                                          int n_rows = -1);
  vil_image_view<double> every_nth_col_of(const vil_image_view<double>& src,
                                          unsigned n, unsigned offset = 0,
                                          int n_cols = -1);

  void reflect_rows(vil_image_view<double>& image, 
                    unsigned top_row, unsigned n_src_rows, 
                    unsigned n_dest_rows);
  void reflect_columns(vil_image_view<double>& image, 
                       unsigned top_row, unsigned n_src_cols, 
                       unsigned n_dest_cols);

  void filter_rows(const vil_image_view<double>& src,
                   vil_image_view<double>& dest, 
                   const vcl_vector<double>& kernel);
  void filter_columns(const vil_image_view<double>& src,
                      vil_image_view<double>& dest, 
                      const vcl_vector<double>& kernel);

  void filt_rows(vil_image_view<double>& src,
                 vil_image_view<double>& dest, 
                 const vcl_vector<double>& kernel);
  void filt_cols(vil_image_view<double>& src,
                 vil_image_view<double>& dest, 
                 const vcl_vector<double>& kernel);

  void rowdfilt(const vil_image_view<double>& src,
                vil_image_view<double>& dest,
                const vcl_vector<double>& kernel);
  void coldfilt(const vil_image_view<double>& src,
                vil_image_view<double>& dest,
                const vcl_vector<double>& kernel);

  void q2c(const vil_image_view<double>& img,
           vil_image_view< vcl_complex<double> >& p_plus_q,
           vil_image_view< vcl_complex<double> >& p_minus_q);

  vcl_complex<double> interpolate(vil_image_view< vcl_complex<double> > z_in,
                                  unsigned row, unsigned column, unsigned level,
                                  vcl_vector<double> weights) const;

  //: Clear all filters
  void clear_filters();
  void reflect_filter(vcl_vector<double>& f);
  void set_dummy_filters();


  //: Set biorthogonal filter
  void set_biort_near_sym_b_bp();
  void set_biort_antonini();

  //: Set q-shift filter
  void set_qshift_qshift_b_bp();

  //
  //  Private variables
  //

  //: The image to be decomposed
  vil_image_view<double> image_;

  //: Number of levels of decomposition
  unsigned n_levels_;
  unsigned level_;

  //: Biorthogonal/q-shift filters
  vcl_vector<double> g0o_, g0a_, g0b_;
  vcl_vector<double> g1o_, g1a_, g1b_;
  vcl_vector<double> g2o_, g2a_, g2b_;
  vcl_vector<double> h0o_, h0a_, h0b_;
  vcl_vector<double> h1o_, h1a_, h1b_;
  vcl_vector<double> h2o_, h2a_, h2b_;

  //: Size of largest filter
  unsigned max_filter_size_;

  //: Workspace images
  vil_image_view<double> im0_;
  vil_image_view<double> work_im_;
  vil_image_view<double> im1_;

  //: Views of workspace images
  vil_image_view<double> src_, dest_;

  unsigned src_left_, src_top_, src_ni_, src_nj_;

  //: Decomposition
  vcl_vector< vil_image_view< vcl_complex<double> > > tree_;

  //: Biorthogonal thingy
  vcl_vector<double> biort_;

  //: QShift thingy
  vcl_vector<double> qshift_;
};

#endif // vil_dtcwt_h_