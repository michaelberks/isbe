#include "ncm_qtool_tmp.h"

#include <QFileDialog>

#include <vcl_cmath.h>

#include <vnl/vnl_random.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/algo/vil_line_filter.h>
#include <vil/algo/vil_suppress_non_max.h>
#include <vil/algo/vil_suppress_non_plateau.h>
#include <vil/algo/vil_threshold.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_corners.h>
#include <vil/algo/vil_dog_filter_5tap.h>
#include <vil/algo/vil_gauss_reduce.h>

#include <nailfold/vil_gaussian_derivatives.h>
#include <nailfold/vil_suppress_non_max_dir.h>
#include <nailfold/vil_find_connected_components.h>

#include <qcore/qcore_convert_image.h>

ncm_qtool_tmp::ncm_qtool_tmp(QWidget *parent, Qt::WFlags flags)
    : QMainWindow(parent, flags)
{
  // setup the UI
  ui.setupUi(this);

  // set contrast and brightness
  brightness_ = 0.0f;
  contrast_ = 0.0f;

  lower_limit_ = 0;
  upper_limit_ = 255;

  // tell graphics view to use this scene
  ui.graphicsView->setScene(&scene_);
  ui.graphicsView->reverse_zoom();

  // Never show scrollbars
  //   The reason for this is that when you switch between the full resolution
  //   image and low resolution 'quick preview' image, the scrollbars interfere
  //   with the coordinate frames. It should be possible to compute the true
  //   correction with scrollbars on, but for now I just switch them off
  ui.graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  ui.graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  // define default colormap
  for (unsigned i = 0; i < 256; i++)
    raw_colour_table_.push_back(qRgb(i,i,i));
  for (unsigned i = 0; i < 256; i++)
    component_colour_table_.push_back(qRgb(0,0,0));
}

ncm_qtool_tmp::~ncm_qtool_tmp()
{
}

void ncm_qtool_tmp::on_actionDebug_activated()
{
}

void ncm_qtool_tmp::update_raw_colour_table()
{
  raw_colour_table_.clear();
  for (int i = 0; i < lower_limit_; i++)
    raw_colour_table_.push_back(qRgb(0,0,0));

  int denom = upper_limit_ - lower_limit_;
  for (int i = lower_limit_; i < upper_limit_; i++)
  {
    int grey = (255 * (i - lower_limit_)) / denom;
    raw_colour_table_.push_back(qRgb(grey,grey,grey));
  }

  for (int i = upper_limit_; i < 256; i++)
    raw_colour_table_.push_back(qRgb(255,255,255));
}

void ncm_qtool_tmp::update_component_colour_table(int value)
{
  int n = component_colour_table_.size();
  for (int i = 0; i <= value; i++)
    component_colour_table_[i] = qRgb(0,0,0);
  for (int i = value+1; i < n; i++)
    component_colour_table_[i] = qRgb(255,255,255);
}

void ncm_qtool_tmp::redraw_raw()
{
  // convert back to QPixmap and display
  QImage qimg;

  //// compute normalized image if necessary
  //if ( (vcl_abs(contrast_-1.0f) > 1e-3) &&
  //     (vcl_abs(brightness_) > 1e-3) )
  //{
  //  normalized_.deep_copy(image_);
  //  vil_math_scale_and_offset_values(normalized_, 1.0+contrast_, brightness_);
  //  vil_image_view<vxl_byte> byte_image;
  //  vil_convert_cast(normalized_, byte_image);
  //  qcore_convert_image(qimg, byte_image);
  //}

  // compute normalized image if necessary
  //if ((lower_limit_>0) || (upper_limit_<255))
  //{
  //  vil_image_view<vxl_byte> byte_image;
  //  vil_convert_stretch_range_limited(image_, byte_image, 
  //                                    static_cast<vxl_byte>(lower_limit_), 
  //                                    static_cast<vxl_byte>(upper_limit_));
  //  qcore_convert_image(qimg, byte_image);
  //}
  //else
    qcore_convert_image(qimg, image_);

  // apply current colour table (colormap)
  qimg.setColorTable(raw_colour_table_);

  // clear scene, add the image and set the scene size accordingly
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
  scene_.setSceneRect(0, 0, qimg.width(), qimg.height());
}

void ncm_qtool_tmp::quick_preview_raw()
{
  // convert back to QPixmap and display
  QImage qimg;

  // compute normalized image if necessary
  if ((lower_limit_>0) || (upper_limit_<255))
  {
    vil_image_view<vxl_byte> byte_image;
    vil_convert_stretch_range_limited(small_image_, byte_image, 
                                      static_cast<vxl_byte>(lower_limit_), 
                                      static_cast<vxl_byte>(upper_limit_));
    qcore_convert_image(qimg, byte_image);
  }
  else
    qcore_convert_image(qimg, small_image_);

  // clear scene, add the image and set the scene size accordingly
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
  scene_.setSceneRect(0, 0, qimg.width(), qimg.height());
}


void ncm_qtool_tmp::on_actionLoad_image_activated()
{
	QString selfilter;
	QString fileName = QFileDialog::getOpenFileName(
    this, "Open file", "u:/projects/nailfold/images",
    tr("All Images (*.jpg *.png *.bmp)" ), // filter
    &selfilter );

  if (fileName.isEmpty())
    return;

  vil_image_view<vxl_byte> loaded_image;
  loaded_image = vil_load(fileName.toStdString().c_str());

  // if image is RGB then convert to greyscale
  vil_image_view<vxl_byte> greyscale_image;
  if (loaded_image.nplanes() > 1)
    vil_convert_planes_to_grey(loaded_image,greyscale_image);
  else
    greyscale_image.deep_copy(loaded_image);

  // Crop top row of pixels
  //   Nailfold images typically come with a row of black pixels along the top
  image_ = vil_crop(greyscale_image, 0, greyscale_image.ni(), 
                                     1, greyscale_image.nj()-1);

  // compute preview image (subsampled)
  vil_image_view<vxl_byte> workim;
  vil_gauss_reduce_121(image_, workim);
  vil_gauss_reduce_121(workim, small_image_);

  // define threshold for 'close to white'
  unsigned const thresh = 256;

  // compute mean of pixels that are not close to white
  double sum = 0;
  unsigned count = 0;
  for (unsigned i = 0; i < image_.ni(); ++i)
    for (unsigned j = 0; j < image_.nj(); ++j)
      if (image_(i,j) <= thresh)
      {
        sum += static_cast<double>(image_(i,j));
        ++count;
      }
  double mean = sum / count;

  // replace close-to-white pixels with mean value
  for (unsigned i = 0; i < image_.ni(); ++i)
    for (unsigned j = 0; j < image_.nj(); ++j)
      if (image_(i,j) > thresh)
        image_(i,j) = mean;

  // convert to a Qt image and display
  QImage qimg;
  qcore_convert_image(qimg, image_);

  qimg.setColorTable(raw_colour_table_);
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
  scene_.setSceneRect(0, 0, qimg.width(), qimg.height());
  ui.graphicsView->fitInView(scene_.sceneRect(),Qt::KeepAspectRatio);
}


void ncm_qtool_tmp::on_actionPan_activated()
{
  ui.graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
}

void ncm_qtool_tmp::on_actionMarkup_activated()
{
  ui.graphicsView->setDragMode(QGraphicsView::NoDrag);
}

void ncm_qtool_tmp::on_actionReset_activated()
{
  redraw_raw();
}

void ncm_qtool_tmp::on_actionSmooth_activated()
{
  vil_gauss_filter_2d(image_, processed_, 2, 8);

  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, processed_);
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_actionFind_lines_activated()
{
  vil_image_view<vxl_byte> line_dir;
  vil_image_view<float> line_str;

  vil_line_filter<vxl_byte> line_filter;
  line_filter.dark_lines_3x3(line_dir, line_str, image_, 0.0001f);

  vil_convert_stretch_range(line_str, line_str, 0, 255);

  vil_image_view<vxl_byte> b_line_str;
  vil_convert_cast(line_str, b_line_str);

  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, b_line_str);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_actionGaussian_deriv_activated()
{
  vil_image_view<double> line_strength;
  vil_image_view<double> orientation;

  // gaussian 2nd derivative processed image
  vil_gaussian_2nd_derivative(image_, 
                              line_strength, orientation);
  
  vil_image_view<vxl_byte> byte_image;
  vil_convert_stretch_range(line_strength, byte_image);
  
  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, byte_image);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_actionMaxima_activated()
{
  vil_image_view<double> line_strength;
  vil_image_view<double> orientation;

  // gaussian 2nd derivative processed image
  vil_gaussian_2nd_derivative(image_, line_strength, orientation,
                              4.0, 24);

  // suppress non-maximal points
  vil_image_view<bool> bool_peaks(image_.ni(), image_.nj());
  vil_suppress_non_max_dir(line_strength, orientation, bool_peaks);

  // label connected components
  vil_find_connected_components(bool_peaks, 
                                component_label_image_, component_sizes_);

  // compute size of biggest component
  unsigned max_size = 0;
  for (unsigned i = 0; i < component_sizes_.size(); i++)
    if (component_sizes_[i] > max_size)
      max_size = component_sizes_[i];

  // Replace label with size of component so that we can remove components by
  // manipulating the colour table rather than the image
  // - Since only 8-bit indexing is supported, we bin the component sizes into
  //   256 bins
  component_size_image_.set_size(component_label_image_.ni(),
                                 component_label_image_.nj());
  for (unsigned j = 0; j < component_label_image_.nj(); j++)
    for (unsigned i = 0; i < component_label_image_.ni(); i++)
      if (component_label_image_(i,j) > 0)
      {
        unsigned real_size = component_sizes_[component_label_image_(i,j)-1];
        component_size_image_(i,j) = (255 * real_size) / max_size;
      }
      else
        component_size_image_(i,j) = 0;

  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, component_size_image_);
  update_component_colour_table(0);
  qimg.setColorTable(component_colour_table_);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_actionHarris_activated()
{
  // compute Harris corners
  vil_image_view<float> grad_i, grad_j;
  vil_sobel_3x3(image_, grad_i, grad_j);
  vil_image_view<float> corners;
  vil_corners(grad_i, grad_j, corners);

  // stretch image range to [0..255]  
  vil_image_view<vxl_byte> byte_image;
  vil_convert_stretch_range(corners, byte_image);
  
  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, byte_image);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_actionDoG_activated()
{
  // convert image from byte to float
  vil_image_view<float> float_image;
  vil_convert_cast(image_,float_image);

  // compute difference of gaussian image
  vil_image_view<float> smooth_im;
  vil_image_view<float> dog_im;
  vil_dog_filter_5tap(float_image, smooth_im, dog_im, 3.0);

  // stretch image range to [0..255]  
  vil_image_view<vxl_byte> byte_image;
  vil_convert_stretch_range(dog_im, byte_image);
  
  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, byte_image);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void ncm_qtool_tmp::on_brightnessSlider_valueChanged(int value)
{
  brightness_ = static_cast<float>(value)*0.01*255;
  redraw_raw();
}

void ncm_qtool_tmp::on_contrastSlider_valueChanged(int value)
{
  contrast_ = static_cast<float>(value)*0.01;
  redraw_raw();
}

void ncm_qtool_tmp::on_lowerLimitSlider_valueChanged(int value)
{
  // upper limit on slider bar value
  int limit = 0;
  if (ui.lockSliders->isChecked())
    limit = lower_limit_ + (255 - upper_limit_);
  else
    limit = upper_limit_;

  // if we hit upper limit then set value to this limit and move the slider to
  // the corresponding position
  if (value > limit)
  {
    value = limit;
    ui.lowerLimitSlider->setValue(value);
  }

  // move upper limit sliderbar if need be
  if (ui.lockSliders->isChecked())
  {
    upper_limit_ = upper_limit_ + (value - lower_limit_);
    lower_limit_ = value;
    if (ui.upperLimitSlider->value() != upper_limit_)
      ui.upperLimitSlider->setValue(upper_limit_);
    else
      return;
  }
  else
    lower_limit_ = value;

  update_raw_colour_table();
  redraw_raw();
}

void ncm_qtool_tmp::on_upperLimitSlider_valueChanged(int value)
{
  // upper limit on slider bar value
  int limit = 0;
  if (ui.lockSliders->isChecked())
    limit = upper_limit_ - (lower_limit_);
  else
    limit = lower_limit_;

  // if we hit upper limit then set value to this limit and move the slider to
  // the corresponding position
  if (value < limit)
  {
    value = limit;
    ui.upperLimitSlider->setValue(value);
  }

  // move lower limit sliderbar if need be
  if (ui.lockSliders->isChecked())
  {
    lower_limit_ = lower_limit_ + (value - upper_limit_);
    upper_limit_ = value;
    if (ui.lowerLimitSlider->value() != lower_limit_)
      ui.lowerLimitSlider->setValue(lower_limit_);
    else
      return;
  }
  else
    upper_limit_ = value;

  update_raw_colour_table();
  redraw_raw();
}

void threshold_sizes(const vil_image_view<int>& labels_image,
                     const vcl_vector<unsigned>& size_vector,
                     vil_image_view<bool>& bool_image,
                     const unsigned size_threshold)
{
  bool_image.set_size(labels_image.ni(), labels_image.nj());

  const int* lbl_row_ptr = labels_image.top_left_ptr();
  bool* out_row_ptr = bool_image.top_left_ptr();

  const unsigned ni = labels_image.ni();
  const unsigned nj = labels_image.nj();

  for (unsigned j = 0; j < nj; j++)
  {
    const int* lbl_ptr = lbl_row_ptr;
    bool* out_ptr = out_row_ptr;

    for (unsigned i = 0; i < ni; i++)
    {
      if (*lbl_ptr > 0)
        *out_ptr = (size_vector[(*lbl_ptr)-1] >= size_threshold);
      else
        *out_ptr = false;

      ++lbl_ptr;
      ++out_ptr;
    }

    lbl_row_ptr += labels_image.jstep();
    out_row_ptr += bool_image.jstep();
  }
}

void ncm_qtool_tmp::on_compSizeSlider_valueChanged(int value)
{
  //vil_image_view<bool> big_components;

  //threshold_sizes(component_label_image_, component_sizes_, 
  //                big_components,
  //                static_cast<unsigned>(value));

  // stretch range and put in vxl_byte image
  vil_image_view<vxl_byte> byte_image;
  //vil_convert_stretch_range(big_components, byte_image);
  //vil_convert_cast(component_sizes, byte_image);

  // convert back to QPixmap and display
  QImage qimg;
  qcore_convert_image(qimg, component_size_image_);
  //qcore_convert_image(qimg, byte_image);

  update_component_colour_table(value);
  qimg.setColorTable(component_colour_table_);

  // add processed image to scene
  scene_.clear();
  scene_.addPixmap(QPixmap::fromImage(qimg));
}

void show_QRect(QRectF r)
{
  vcl_cout << "QRectF: "
           << "(" << r.left() << "," << r.top() << ")"
           << " -> "
           << "(" << r.right() << "," << r.bottom() << ")"
           << vcl_endl;
}

void show_QPoint(QPointF p)
{
  vcl_cout << "QPointF: "
           << "(" << p.x() << "," << p.y() << ")"
           << vcl_endl;
}

void ncm_qtool_tmp::on_pushButton_clicked()
{
  QPoint view_centre = ui.graphicsView->childrenRect().center();

  preview_scene_centre_ = ui.graphicsView->mapToScene(view_centre);
  show_QPoint(preview_scene_centre_);

  ui.graphicsView->centerOn(preview_scene_centre_);

  show_QPoint(ui.graphicsView->mapToScene(
    ui.graphicsView->childrenRect().center()) );

  vcl_cout << vcl_endl;
}

//
//
void ncm_qtool_tmp::on_lowerLimitSlider_sliderPressed()
{
  //// translate to centre of low-res preview image
  //QPoint view_centre = ui.graphicsView->contentsRect().center();
  //preview_scene_centre_ = ui.graphicsView->mapToScene(view_centre);
  //preview_scene_centre_ += QPointF(2.0, 2.0);
  //ui.graphicsView->centerOn((preview_scene_centre_/4.0));

  //// zoom in to account for scale difference
  //ui.graphicsView->scale(4.0, 4.0);

  //// draw low-res preview image
  //quick_preview_raw();
}

void ncm_qtool_tmp::on_upperLimitSlider_sliderPressed()
{
  //// translate to centre of low-res preview image
  //QPoint view_centre = ui.graphicsView->contentsRect().center();
  //preview_scene_centre_ = ui.graphicsView->mapToScene(view_centre);
  //preview_scene_centre_ += QPointF(2.0, 2.0);
  //ui.graphicsView->centerOn((preview_scene_centre_/4.0));

  //// zoom in to account for scale difference
  //ui.graphicsView->scale(4.0, 4.0);

  //// draw low-res preview image
  //quick_preview_raw();
}

void ncm_qtool_tmp::on_lowerLimitSlider_sliderReleased()
{
  //redraw_raw();
  //preview_scene_centre_ -= QPointF(2.0, 2.0);
  //ui.graphicsView->centerOn(preview_scene_centre_);
  //ui.graphicsView->scale(0.25, 0.25);
}

void ncm_qtool_tmp::on_upperLimitSlider_sliderReleased()
{
  //redraw_raw();
  //preview_scene_centre_ -= QPointF(2.0, 2.0);
  //ui.graphicsView->centerOn(preview_scene_centre_);
  //ui.graphicsView->scale(0.25, 0.25);
}
