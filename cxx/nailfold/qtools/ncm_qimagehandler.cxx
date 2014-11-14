//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qimagehandler.h"

#include "ncm_qglobal.h"

#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>
#include <vil/vil_fill.h>

#include <vil/algo/vil_threshold.h>
#include <vil/algo/vil_structuring_element.h>
#include <vil/algo/vil_binary_dilate.h>

#include <vul/vul_file.h>

#include <mbl/mbl_index_sort.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/vil_gaussian_derivatives.h>
#include <nailfold/vil_suppress_non_max_dir.h>
//#include <nailfold/vil_find_connected_components.h>

#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>

//
//  Public methods
//

ncm_qimagehandler::ncm_qimagehandler(QWidget* parent)
: QObject(parent), // base class construction
  raw_image_(QImage()),
  line_image_(QImage()),
  edge_image_(QImage()),
  raw_pixmap_(QPixmap()),
  line_pixmap_(QPixmap()),
  edge_pixmap_(QPixmap()),
  min_contrast_(0),
  max_contrast_(255),
  line_size_threshold_(1),
  invalid_(0),
  image_type_(TypeUndefined)
{
  ++invalid_;

  // initialize colour tables
  update_raw_colour_table();
  update_centreline_colour_table();
  update_edge_colour_table();

  // use binned size image
  conn_comp_.set_binning_strategy(vil_connected_components::BinningEqualCounts);

  --invalid_;
}

//
//: Destructor
ncm_qimagehandler::~ncm_qimagehandler()
{
}

void ncm_qimagehandler::clear()
{
  raw_vxl_image_.clear();

  raw_image_ = line_image_ = edge_image_ = QImage();
  raw_pixmap_ = line_pixmap_ = edge_pixmap_ = QPixmap();

  conn_comp_.clear();

  line_pixels_.clear();
  edge_pixels_.clear();
}

bool ncm_qimagehandler::isValid() const
{
  return (invalid_ == 0);
}

int ncm_qimagehandler::min_contrast() const
{
  return min_contrast_;
}
int ncm_qimagehandler::max_contrast() const
{
  return max_contrast_;
}

bool ncm_qimagehandler::image_is_large() const
{
  return (raw_image_.width() > 1000);
}

void ncm_qimagehandler::set_contrast(int min_contrast, int max_contrast,
                                     bool auto_update /* = true */)
{
  if (min_contrast == min_contrast_ && 
      max_contrast == max_contrast_)
    return;

  ++invalid_;

  if (min_contrast < 0)
    min_contrast_ = 0;
  else if (min_contrast > 255)
    min_contrast_ = 255;
  else
    min_contrast_ = min_contrast;

  if (max_contrast < 0)
    max_contrast_ = 0;
  else if (max_contrast > 255)
    max_contrast_ = 255;
  else
    max_contrast_ = max_contrast;

  if (auto_update)
    update_raw_colour_table();

  --invalid_;
}
void ncm_qimagehandler::set_min_contrast(int min_contrast,
                                  bool auto_update /* = true */)
{
  if (min_contrast == min_contrast_)
    return;

  ++invalid_;

  if (min_contrast < 0)
    min_contrast = 0;
  else if (min_contrast > 255)
    min_contrast = 255;

  min_contrast_ = min_contrast;
  if (auto_update)
    update_raw_colour_table();

  --invalid_;
}
void ncm_qimagehandler::set_max_contrast(int max_contrast,
                                  bool auto_update /* = true */)
{
  if (max_contrast == max_contrast_)
    return;

  ++invalid_;

  if (max_contrast < 0)
    max_contrast = 0;
  else if (max_contrast > 255)
    max_contrast = 255;

  max_contrast_ = max_contrast;
  if (auto_update)
    update_raw_colour_table();

  --invalid_;
}

void ncm_qimagehandler::set_contrast_from(QRect image_rect)
{
  // Do nothing if this is a (colour) dermatogram.
  if (image_type_ == TypeDermatogram)
    return;

  ++invalid_;

  const unsigned i0 = vcl_max<int>(image_rect.left(), 0);
  const unsigned ni = vcl_min<int>(image_rect.width(), raw_vxl_image_.ni()-i0);
  const unsigned j0 = vcl_max<int>(image_rect.top(), 0);
  const unsigned nj = vcl_min<int>(image_rect.height(), raw_vxl_image_.nj()-j0);

  if (ni == 0 || nj == 0)
  {
    --invalid_;
    return;
  }

  vil_image_view<vxl_byte> view = vil_crop(raw_vxl_image_, i0, ni, j0, nj);
      
  vxl_byte min_value = 255;
  vxl_byte max_value = 0;

  // Customize vil_math_value_range() here to ignore the brightest value
  // pixels that cause sudden jumps in image brightness when using
  // autocontrast.
  unsigned n_pixels = 0;
  for (unsigned j=0; j<nj; ++j)
  {
    for (unsigned i=0; i<ni; ++i)
    {
      const vxl_byte pixel = view(i,j);

      if (pixel == 255 || pixel == 0)
        continue;
      if (pixel < min_value)
        min_value = pixel;
      if (pixel > max_value)
        max_value = pixel;
      ++n_pixels;
    }
  }

	if (min_value == 255) //Image contains no values between 0 and 255
		min_value = 0;

	if (max_value == 0) //Image contains no values between 0 and 255
		max_value = 255;
  
  set_contrast(min_value, max_value);

  --invalid_;
}

void ncm_qimagehandler::set_line_size_threshold(int line_size_threshold)
{
  ++invalid_;

  // range checking
  if (line_size_threshold < 1)
    line_size_threshold_ = 1;
  else if (line_size_threshold > 255)
    line_size_threshold_ = 255;
  else
    line_size_threshold_ = line_size_threshold;

  update_centreline_colour_table();

  --invalid_;
}

//
//
void ncm_qimagehandler::set_filename(const vcl_string& filename)
{
  filename_ = filename;
}
vcl_string ncm_qimagehandler::filename() const
{
  return filename_;
}

unsigned ncm_qimagehandler::width() const
{
  if (raw_image_.isNull())
    return 0;
  else
    return raw_image_.width();
}
unsigned ncm_qimagehandler::height() const
{
  return raw_image_.height();
}

QPixmap const& ncm_qimagehandler::raw_pixmap() const
{
  return raw_pixmap_;
}
QPixmap const& ncm_qimagehandler::line_pixmap() const
{
  return line_pixmap_;
}
QPixmap const& ncm_qimagehandler::edge_pixmap() const
{
  return edge_pixmap_;
}

void ncm_qimagehandler::load_image()
{
  ++invalid_;
  load_image_from(filename_);
  --invalid_;
}

//
//: Load an image from a given filename
void ncm_qimagehandler::load_image_from(const vcl_string& filename)
{
  ++invalid_;

  // load the image
  vil_image_view<vxl_byte> image;
  set_status("Loading image...");
  image = vil_load(filename.c_str());

  // Determine the image type
  vcl_string image_basename = vul_file::basename(filename);
  image_basename = vul_file::strip_extension(image_basename);
  if (image_basename[image_basename.length()-1] == 'd' &&
      image.nplanes() == 3)
    set_image_type(TypeDermatogram);
  else
    set_image_type(TypeCapillarogram);

  set_status("");
  set_image(image);

  --invalid_;
}

//
//: Set nailfold image
void ncm_qimagehandler::set_image(const vil_image_view<vxl_byte>& image)
{
  ++invalid_;

  switch (image_type_)
  {
    case TypeCapillarogram:
    {
      // if image is RGB, convert to greyscale
      set_status("Flattening image...", 0);
      vil_image_view<vxl_byte> greyscale_image;
      if (image.nplanes() > 1)
        //vil_convert_planes_to_grey(image, greyscale_image);
        greyscale_image = vil_plane(image, 1); // green channel
      else
        greyscale_image = image;

      // remove top row of pixels
      //   Nailfold images typically come with a row of black pixels along the top
      //vil_image_view<vxl_byte> raw_vxl_image;
      raw_vxl_image_ = vil_crop(greyscale_image, 0, greyscale_image.ni(), 
                                                 1, greyscale_image.nj()-1);

      // Create a mask of the background pixels.
      vil_image_view<bool> bg_mask(raw_vxl_image_.ni(), raw_vxl_image_.nj());
      bg_mask.fill(false);
      for (unsigned i = 0; i < bg_mask.ni(); ++i)
        for (unsigned j = 0; j < bg_mask.nj(); ++j)
          if (raw_vxl_image_(i,j) == 255)
          {
            bg_mask(i,j) = true;
          }

      // Dilate the background mask to include antialiased border pixels.
      vil_image_view<bool> bg_mask_dilated(bg_mask.ni(), bg_mask.nj());
      bg_mask_dilated.fill(false);

      vil_structuring_element strel;
      strel.set_to_disk(3);

      vil_binary_dilate(bg_mask, bg_mask_dilated, strel);

      // Fill the region of the image covered by the mask (i.e. background and
      // border pixels) with pure white (255).
      vil_fill_mask<vxl_byte>(raw_vxl_image_, bg_mask_dilated, 255);

      // convert to QImage and set scene elements
      qcore_convert_image(raw_image_, raw_vxl_image_);

      break;
    }

    case TypeDermatogram:
    {
      qcore_convert_image(raw_image_, image);
      break;
    }

    default:
      break;
  }

  // TODO: 1st and 2nd derivatives are different in that one gives the
  // orientation normal to the edge whereas the other gives the orientation
  // along the line. This should be corrected and a normals() function
  // added to the class (when I've converted it to a class, that is).

  // Apply colourmaps
  apply_raw_colour_table();

  set_status("", -1.0); // reset

  emit imageLoaded();

  --invalid_;
}

//
//: Compute Gaussian second derivatives
void ncm_qimagehandler::find_ridges()
{
  ++invalid_;

  if (line_image_.width() > 0)
  {
    // already done
    --invalid_;
    return;
  }

  // compute Gaussian 2nd derivatives
  vil_image_view<double> g2d_strength, g2d_orientation;
  set_status("Filtering...", 0.2);
  vil_gaussian_2nd_derivative(raw_vxl_image_, 
                              g2d_strength, g2d_orientation,
                              4.0);

  // create workspace images for peaks and their maxima
  vil_image_view<double> double_image(raw_vxl_image_.ni(), raw_vxl_image_.nj());
  vil_image_view<bool> bool_image(raw_vxl_image_.ni(), raw_vxl_image_.nj());

  // suppress non-maximal points
  vil_image_view<double>& peaks_g2d = double_image;
  set_status("Finding maxima...", 0.6);
  vil_suppress_non_max_dir(g2d_strength, g2d_orientation,
                           peaks_g2d, false);

  // threshold to get only those with a response > threshold
  vil_image_view<bool>& bool_peaks_g2d = bool_image;
  set_status("Thresholding...", 0.75);
  vil_threshold_above(peaks_g2d, bool_peaks_g2d, 1e-6);

  // compute connected components of gaussian 2nd derivative lines
  set_status("Finding connected components...", 0.8);
  conn_comp_.set_image(bool_peaks_g2d);

  // convert to QImage and set scene elements
  set_status("Converting...", 0.9);
  qcore_convert_image(line_image_, conn_comp_.size_image());

  // Set colourmap and find ridge pixels
  apply_centreline_colour_table();
  get_line_pixels();

  set_status("", -1.0); // reset

  --invalid_;
}

void ncm_qimagehandler::find_edges()
{
  ++invalid_;

  if (edge_image_.width() > 0)
  {
    // already done
    --invalid_;
    return;
  }

  // compute Gaussian 2nd derivatives
  vil_image_view<double> g1d_strength, g1d_orientation;
  set_status("Filtering...", 0.2);
  vil_gaussian_1st_derivative(raw_vxl_image_, 
                              g1d_strength, g1d_orientation,
                              4.0);

  // create workspace images for peaks and their maxima
  vil_image_view<double> double_image(raw_vxl_image_.ni(), raw_vxl_image_.nj());
  vil_image_view<bool> bool_image(raw_vxl_image_.ni(), raw_vxl_image_.nj());

  // suppress non-maximal points
  vil_image_view<double>& peaks_g1d = double_image;
  set_status("Finding maxima...", 0.6);
  vil_suppress_non_max_dir(g1d_strength, g1d_orientation,
                           peaks_g1d, true);

  // threshold to get only those with a response > threshold
  vil_image_view<bool>& bool_peaks_g1d = bool_image;
  set_status("Thresholding...", 0.75);
  vil_threshold_above(peaks_g1d, bool_peaks_g1d, 1e-6);

  // compute connected components of gaussian 2nd derivative lines
  set_status("Finding connected components...", 0.8);
  conn_comp_.set_image(bool_peaks_g1d);

  // convert to QImage and set scene elements
  set_status("Converting...", 0.9);
  qcore_convert_image(edge_image_, conn_comp_.size_image());

  // Set colourmap and find edge pixels
  apply_edge_colour_table();
  get_edge_pixels();

  --invalid_;
}

//
//: Find pixels that are labelled as a detected line and are above the selected
//  threshold
void ncm_qimagehandler::get_line_pixels()
{
  ++invalid_;
  get_pixmap_pixels(line_pixmap_, line_pixels_);
  --invalid_;
}

//
//: Find pixels that are labelled as a detected edge and are above the selected
//  threshold
void ncm_qimagehandler::get_edge_pixels()
{
  ++invalid_;
  get_pixmap_pixels(edge_pixmap_, edge_pixels_);
  --invalid_;
}

ncm_qimagehandler::imageType ncm_qimagehandler::image_type() const
{
  return image_type_;
}

//
//: Find nearest visible line pixel to the given (x,y) coordinate
QPoint ncm_qimagehandler::line_pixel_nearest_to(const QPointF& pos) const
{
  return pixel_nearest_to(pos, line_pixels_);
}

QPoint ncm_qimagehandler::edge_pixel_nearest_to(const QPointF& pos) const
{
  return pixel_nearest_to(pos, edge_pixels_);
}

//
//: Change colormap for raw images to effect a change in contrast
void ncm_qimagehandler::update_raw_colour_table()
{
  ++invalid_;

  raw_colour_table_.resize(256);

  // fill everything from first to min_contrast_ with 0 (black)
  for (int i = 0; i < min_contrast_; i++)
    raw_colour_table_[i] = NcmQt::black;

  // fill everything from max_contrast_ to end with 255 (white)
  for (int i = max_contrast_; i < 256; i++)
    raw_colour_table_[i] = NcmQt::white;

  // fill values in between with a linear gradient
  int denom = max_contrast_ - min_contrast_;
  for (int i = min_contrast_; i < max_contrast_; i++)
  {
    int intensity = (255 * (i - min_contrast_)) / denom;
    raw_colour_table_[i] = qRgb(intensity,intensity,intensity);
  }

  apply_raw_colour_table();

  --invalid_;
}

void ncm_qimagehandler::apply_raw_colour_table()
{
  ++invalid_;

  if (!raw_image_.isNull())
  {
    switch (image_type_)
    {
      case TypeCapillarogram:
        // update QImage and its associated pixmap
        raw_image_.setColorTable(raw_colour_table_);
        raw_pixmap_ = QPixmap::fromImage(raw_image_);
        break;

      case TypeDermatogram:
      {
        raw_pixmap_ = QPixmap::fromImage(raw_image_);
        break;
      }

      default:
        break;
    }

    emit rawImageChanged();
  }

  --invalid_;
}

//
//: Change colormap for vessel centrelines to hide those below the threshold
void ncm_qimagehandler::update_centreline_colour_table()
{
  ++invalid_;

  centreline_colour_table_.resize(256);

  // set everything below the threshold to 0 (transparent black)
  for (int i = 0; i < line_size_threshold_; i++)
    centreline_colour_table_[i] = NcmQt::transparent;

  // set everything above the threshold to 255 (opaque purple)
  for (int i = line_size_threshold_; i < 256; i++)
    centreline_colour_table_[i] = NcmQt::purple;

  apply_centreline_colour_table();

  --invalid_;
}

void ncm_qimagehandler::apply_centreline_colour_table()
{
  ++invalid_;

  if (!line_image_.isNull())
  {
    // update QImage and its associated pixmap
    line_image_.setColorTable(centreline_colour_table_);
    line_pixmap_ = QPixmap::fromImage(line_image_);

    emit lineImageChanged();
  }

  --invalid_;
}

//
//: Change colormap for vessel edges to hide those below the threshold
void ncm_qimagehandler::update_edge_colour_table()
{
  ++invalid_;

  edge_colour_table_.resize(256);

  // set everything below the threshold to 0 (transparent black)
  for (int i = 0; i < line_size_threshold_; i++)
    edge_colour_table_[i] = NcmQt::transparent;

  // set everything above the threshold to 255 (opaque purple)
  for (int i = line_size_threshold_; i < 256; i++)
    edge_colour_table_[i] = NcmQt::yellow;

  apply_edge_colour_table();

  --invalid_;
}

void ncm_qimagehandler::apply_edge_colour_table()
{
  ++invalid_;

  if (!edge_image_.isNull())
  {
    // update QImage and its associated pixmap
    edge_image_.setColorTable(edge_colour_table_);
    edge_pixmap_ = QPixmap::fromImage(edge_image_);

    emit edgeImageChanged();
  }

  --invalid_;
}

void ncm_qimagehandler::set_image_type(ncm_qimagehandler::imageType image_type)
{
  if (image_type_ != image_type)
    emit imageTypeChanged();

  image_type_ = image_type;
}

//
//: Find pixels that are labelled as a detected line and are above the selected
//  threshold
void ncm_qimagehandler::get_pixmap_pixels(const QPixmap& pixmap, 
                                   vcl_vector<QPoint>& pixel_list)
{
  ++invalid_;

  // clear current list
  pixel_list.resize(0);

  // if no pixmap computed then do nothing (return)
  if (pixmap.isNull())
  {
    --invalid_;
    return;
  }

  QImage thresholded_image = pixmap.toImage();
  for (int y = 0; y < thresholded_image.height(); ++y)
    for (int x = 0; x < thresholded_image.width(); ++x)
    {
      if (qAlpha(thresholded_image.pixel(x,y)) > 0)
        pixel_list.push_back(QPoint(x,y));
    }

  --invalid_;
}

//
//: Find nearest visible pixel to the given (x,y) coordinate
QPoint ncm_qimagehandler::pixel_nearest_to(const QPointF& pos,
                                    const vcl_vector<QPoint>& pixel_list) const
{
  // search through list for nearest point
  // could be done more effectively by storing nonzero pixels as a
  // vector< vector<QPoint> > with one vector per row (or column) of the image
  // the list could then be searched starting from the current row/col,
  // propagating the search raster outward

  QPoint nearest_point(0,0);
  int min_distance = -1;

  const unsigned n_pixels = pixel_list.size();
  for (unsigned i = 0; i < n_pixels; ++i)
  {
    int d = (pos - pixel_list[i]).manhattanLength();

    if ((min_distance == -1) || (d < min_distance))
    {
      nearest_point = pixel_list[i];
      min_distance = d;
    }
  }

  return nearest_point;
}

//
//  Private methods
//

//
//: Set status string
void ncm_qimagehandler::set_status(const vcl_string& status, 
                            double progress /* = 0.0 */)
{
  status_ = status;
  emit statusChanged(QString(status_.c_str()), progress);
}
