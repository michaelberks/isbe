#include "ncm_qstitch_main.h"

#include <vcl_iomanip.h>
#include <vcl_algorithm.h>
#include <vcl_sstream.h>
#include <vcl_iostream.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_crop.h>

#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_threshold.h>
#include <vil/algo/vil_gauss_reduce.h>

#include <vul/vul_file_iterator.h>
#include <vul/vul_string.h>
#include <vul/vul_file.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/flow/ncm_flow_field.h>
#include <nailfold/flow/ncm_flow_estimator_cholmod.h>

#include <QImage>
#include <QPixmap>
#include <QGraphicsPixmapItem>
#include <QFileDialog>

//
//  Global functions
//

//: Find the longest string of consecutive numbers from a list.
vcl_vector<unsigned> find_consecutive(vcl_vector<unsigned> indices)
{
  if (indices.empty())
    return indices;

  vcl_sort(indices.begin(), indices.end());

  unsigned longest_string = 1;
  unsigned string_length = 1;
  unsigned first_index = 0, last_index = 0;
  unsigned best_first = 0, best_last = 0;

  while (first_index + longest_string < indices.size())
  {
    ++last_index;

    if (indices[last_index] == indices[last_index-1]+1)
    {
      ++string_length;

      if (string_length > longest_string)
      {
        longest_string = string_length;
        best_first = first_index;
        best_last = last_index;
      }
    }
    else
    {
      first_index = last_index;
      string_length = 1;
    }
  }

  // Return indices[first_index..last_index]
  indices.resize(best_last+1);
  indices.erase(indices.begin(), indices.begin()+best_first);
  return indices;
}

//
//  Public member functions
//
 
//: Constructor
ncm_qstitch_main::ncm_qstitch_main(QWidget *parent, Qt::WFlags flags)
: border_item_(NULL),
  pixelCountItem_(NULL),
  flowRectItem_(NULL),
  maskItem_(NULL),
  stitcher_(NULL)
{
  // setup the UI
  ui_.setupUi(this);
  ui_.mosaicView->setScene(&scene_);

  initializeStatusbar();
  initializeAlignerThread();

  QObject::connect( &animateTimer_, SIGNAL(timeout()),
                    this, SLOT(onAnimateTimer()) );

  // Hide the debug button in Release code.
  #if NDEBUG
    ui_.debug->hide();
  #endif
}
 
//: Destructor
ncm_qstitch_main::~ncm_qstitch_main()
{
  // Clean up.
  reset();

  // This is the only pointer not automatically handled by a parent object.
  delete stitcher_;
}

 
//
//  Public slots
//

//
//  Protected methods
//

//
//  Private slots
//
 
//: Select the input image folder.
void ncm_qstitch_main::on_actionSelectFolder_triggered()
{
  select_folder();
}
 
//: Select the input image folder.
void ncm_qstitch_main::on_selectFolder_clicked()
{
  select_folder();
}
 
//: Limit the number of frames aligned and stitched into a mosiac.
void ncm_qstitch_main::on_nFramesSpin_valueChanged(int nFrames)
{
  aligner_.setNumToAlign(nFrames);
  ui_.frameSlider->setMaximum(nFrames);
}
 
//: Erase existing displacements in the aligner class.
void ncm_qstitch_main::on_clearDisplacements_clicked()
{
  aligner_.resetDisplacements();
}
 
//: Create a mask so that dirt on the lens (which is in the same place at every
//  frame and can therefore be detected by averaging all the frames) is masked
//  out before alignment.
void ncm_qstitch_main::on_createMask_clicked()
{
  vil_image_view<vxl_byte> img;
  vil_image_view<double> mean_image;

  // Compute the mean image (where dirt should be most obvious).
  for (unsigned i = 0; i < n_images(); ++i)
  {
    img = vil_load(filenames_[i].c_str());
  
    if (mean_image.ni() == 0)
      vil_convert_cast(img, mean_image);
    else
    {
      vil_image_view<double> double_image;
      vil_convert_cast(img, double_image);
      vil_math_image_sum(mean_image, double_image, mean_image);
    }
  }
  vil_math_scale_values(mean_image, 1.0/n_images());

  // Find pixels with high gradient magnitude.
  vil_image_view<double> grad_ij;
  vil_sobel_3x3(mean_image, grad_ij);

  vil_image_view<double> grad_magnitude;
  vil_math_rss(grad_ij, grad_magnitude);

  // Threshold (at a multiple of the mean value) to give a binary mask of 
  // dirt edge pixels.
  double grad_threshold;
  vil_math_mean(grad_threshold, grad_magnitude, 0);

  vil_image_view<bool> mask_image;
  vil_threshold_above(grad_magnitude, mask_image, 2.0*grad_threshold);

  aligner_.set_mask(mask_image);
}
 
//: Show/hide the mask in the scene view.
void ncm_qstitch_main::on_showMask_toggled(bool checked)
{
  if (checked && (aligner_.mask().ni() > 0))
  {
    // Convert from binary to greyscale so that it can be converted to a 
    // pixmap.
    vil_image_view<vxl_byte> grey_image;
    vil_convert_cast(aligner_.mask(), grey_image);

    vil_convert_stretch_range(grey_image, grey_image);
    maskItem_ = scene_.addPixmap(pixmap_from(grey_image));
  }
  else
  {
		//*** this breaks if a scene is present!
    scene_.removeItem(maskItem_);
    delete maskItem_;
    maskItem_ = NULL;
  }
}
 
//: Copy the displacements from the motor positions to the aligner as an
//  initialization. In practice, 'good' sequences don't need this, though
//  sequences with large regions of uniform texture may benefit from it.
void ncm_qstitch_main::on_useAllMotorPositions_clicked()
{
  aligner_.resetDisplacements();
  for (unsigned i = 1; i < n_images(); ++i)
  {
    // Compute displacements from motor positions.
    int di = motor_positions_[i].x() - motor_positions_[i-1].x();
    int dj = motor_positions_[i].y() - motor_positions_[i-1].y();

    aligner_.setDisplacement(i, QPoint(di, dj));
  }
}
 
//: Change the maximum distance allowed between consecutive frames.
void ncm_qstitch_main::on_radiusSpin_valueChanged(int value)
{
  aligner_.set_di_radius(value);
  aligner_.set_dj_radius(value);
}
 
//: Align the frames in the sequence using the fast Sobel method on heavily
//  downsampled images.
void ncm_qstitch_main::on_coarseAlign_clicked()
{
  if (ui_.coarseAlign->isChecked())
  {
		if (aligner_.numToAlign() < 2)
			return;

    aligner_.set_filter(ncm_frame_aligner::filter_sobel);
    aligner_.set_n_levels(2);
    aligner_.set_di_radius(ui_.radiusSpin->value());
    aligner_.set_dj_radius(ui_.radiusSpin->value());

    ui_.coarseAlign->setText("Stop Align");
    start_alignment();
  }
  else
  {
    aligner_thread_.exit();
    ui_.coarseAlign->setText("Start Align");
  }
}
 
//: Align the frames in the sequence using the slower Gaussian 1st derivatives
//  method on non-downsampled images.
void ncm_qstitch_main::on_fineAlign_clicked()
{
  if (ui_.fineAlign->isChecked())
  {
    aligner_.set_filter(ncm_frame_aligner::filter_g1d);
    aligner_.set_n_levels(0);
    ui_.fineAlign->setText("Stop Align");
    start_alignment();
  }
  else
  {
    aligner_thread_.exit();
    ui_.fineAlign->setText("Start Align");
  }
}
 
//: React to two frames having been aligned by adding the more recent image
//  to the scene at its aligned position.
void ncm_qstitch_main::onFramesAligned(
  int di, int dj, 
  double scale, double offset, 
  double mse)
{
  unsigned n_aligned = aligner_.numAligned();

  // Add the pixmap at the origin in the first instance.
  QPixmap pixmap = pixmap_from(aligner_.destination_image(), 1.0, scale, offset);
  frame_items_[n_aligned-1] = add_pixmap(pixmap, QPoint(0,0));

  // Update the positions of all aligned frames.
  update_frameitem_positions(n_aligned-1);

  // Update the status bar.
  updateStatus("Aligning images", 
               static_cast<double>(n_aligned) / n_images());

  // Ask for next frames to be aligned.
  emit readyToAlign();
}
 
//: React to all frames having been aligned.
void ncm_qstitch_main::onAlignmentFinished()
{
  // Update the status bar.
  double dtime = double(vcl_clock() - start_time) / CLOCKS_PER_SEC;
  vcl_stringstream ss;
  ss << "Alignment finished in " << dtime << "s"
     << "(" << n_images() / dtime << " fps)";

  updateStatus(ss.str().c_str(), 0);

  // Stop the alignment thread.
  aligner_thread_.exit();

  // Restore state of align buttons.
  ui_.coarseAlign->setText("Start Align");
  ui_.coarseAlign->setChecked(false);
  ui_.fineAlign->setText("Start Align");
  ui_.fineAlign->setChecked(false);

  // Enable the next control in the sequence
  if (!ui_.fineBox->isEnabled())
    ui_.fineBox->setEnabled(true);
  if (!ui_.blendBox->isEnabled())
    ui_.blendBox->setEnabled(true);

  // Add a border (that will be around the currently 'selected' frame) around 
  // the first frame.
  if (border_item_ == NULL)
  {
    border_item_ = scene_.addRect(frame_items_[0]->boundingRect(),
                                  QPen(QColor(255,0,0)), QBrush());
  }
}
 
//: Change the horizontal displacement of the currently selected frame (which
//  also affects all subsequent frames).
void ncm_qstitch_main::on_horizDispSpin_valueChanged(int dx)
{
  const unsigned frame = ui_.frameSlider->value()-1;

  QPoint d = aligner_.displacement(frame);
  d.setX(dx);
  aligner_.setDisplacement(frame, d);

  // Update subsequent frame positions.
  update_frameitem_positions(/* first = */ frame);
}
 
//: Change the vertical displacement of the currently selected frame (which
//  also affects all subsequent frames).
void ncm_qstitch_main::on_vertDispSpin_valueChanged(int dy)
{
  const unsigned frame = ui_.frameSlider->value()-1;

  QPoint d = aligner_.displacement(frame);
  d.setY(dy);
  aligner_.setDisplacement(frame, d);

  // Update subsequent frame positions.
  update_frameitem_positions(/* first = */ frame);
}
 
//: Copy the displacement from the motor positions for this frame only.
void ncm_qstitch_main::on_useMotorPositions_clicked()
{
  const unsigned frame = ui_.frameSlider->value()-1;

  int di = motor_positions_[frame].x() - motor_positions_[frame-1].x();
  int dj = motor_positions_[frame].y() - motor_positions_[frame-1].y();
  aligner_.setDisplacement(frame, QPoint(di, dj));

  // Update subsequent frame positions.
  update_frameitem_positions(/* first = */ frame);
}
 
//: Change the currently selected frame in response to sliding the slider.
void ncm_qstitch_main::on_frameSlider_valueChanged(int frame)
{
  //if (frame < 2) 
  //  return;

  //vil_image_view<vxl_byte> img;

  //double scale = 0.5;

  //// Vectors are zero-based
  //img = vil_load(filenames_[frame-2].c_str());
  //ui_.previousFrame->setPixmap( pixmap_from(img, scale) );

  //img = vil_load(filenames_[frame-1].c_str());
  //ui_.currentFrame->setPixmap( pixmap_from(img, scale) );

  //if (border_item_ != NULL)
  //  border_item_->setPos(frame_items_[frame-1]->pos());

  //QPoint d = aligner_.displacement(frame-1);
  //ui_.horizDispSpin->setValue(d.x());
  //ui_.vertDispSpin->setValue(d.y());
  //
  //// Label is one-based (more intuitive to humans)
  //ui_.frameTitle->setText("Frame " + QString::number(frame));
}
 
//: Create a mosaic of the aligned images by blending.
void ncm_qstitch_main::on_blendButton_clicked()
{
  // Start the stitcher thread.
  n_stitched_ = 0;
  initializeStitcherThread();
  stitcher_thread_.start();

  // Begin stitching but only for the chosen number of files.
  // (This is not done elegantly but works ok.)
  vcl_vector<vcl_string> filenames = filenames_;
  filenames.resize(n_images());
  vcl_vector<QPoint> displacements = aligner_.displacements();
  displacements.resize(n_images());
  emit makeMosaic(filenames, displacements);
}
 
//: React to the mosaic being updated.
void ncm_qstitch_main::onMosaicUpdated()
{
  // Update the status bar.
  ++n_stitched_;
  updateStatus("Creating mosaic", 
               static_cast<double>(n_stitched_) / n_images());
}
 
//: React to the mosaic being finished by displaying in the scene view.
void ncm_qstitch_main::onMosaicFinished()
{
  // Update the status bar first.
  updateStatus("Mosaic finished", 0);
  stitcher_thread_.exit();

  // Clear the scene and add the mosaic.
  show_as_mosaic();

  // Enable the 'save' button.
  if (!ui_.outputBox->isEnabled())
    ui_.outputBox->setEnabled(true);

  // Begin the timer so that cropped videos are displayed as the mouse moves
  // over the image.
  const double clicks_per_second = 30;
  animateTimer_.setInterval(1000.0 / clicks_per_second);
  animateTimer_.start();
}

//: Show the pixel count image (i.e. a monochrome in which pixels visible in
//  the most number of frames are shown in white, whereas those visible in the
//  least are shown in black).
void ncm_qstitch_main::on_showPixelCount_toggled(bool checked)
{
  if (checked)
  {
    // Create the pixmapitem if need be.
    if (pixelCountItem_ == NULL)
    {
      // Convert from unsigned to double and stretch to range [0,255].
      vil_image_view<double> double_count;
      vil_convert_cast(stitcher_->count_image(), double_count);
      vil_image_view<vxl_byte> img;
      vil_convert_stretch_range(double_count, img);

      // Add a pixmap based on the resulting image.
      QPixmap pixmap = pixmap_from(img);
      pixelCountItem_ = scene_.addPixmap(pixmap);

      // Resize the scene to contain just the mosaic.
      ui_.mosaicView->fitInView(scene_.sceneRect(), Qt::KeepAspectRatio);
    }

    pixelCountItem_->show();
  }
  else
  {
    pixelCountItem_->hide();
  }
}

//: Save the mosaic to a desired filename.
void ncm_qstitch_main::on_saveImage_clicked()
{
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save File"),
    root_dir_.c_str(),
    tr("Images (*.png *.jpg)"));

  if (!fileName.isNull())
    vil_save(stitcher_->mosaic(), fileName.toStdString().c_str());
}
 
//: React to the timer timeout by updating the frame in the respective views.
void ncm_qstitch_main::onAnimateTimer()
{
  // Update the frame in the live view (i.e. showing the region directly
  // below the current mouse position).
  static unsigned vframe_index = 0;

  if (visible_frames_.empty())
    ui_.liveView->setPixmap(QPixmap());
  else
  {
    // Go to the next visible frame.
    if ( (visible_frames_.front() <= vframe_index) &&
                                    (vframe_index < visible_frames_.back()) )
      ++vframe_index;
    else
      vframe_index = visible_frames_.front();

    // Load the frame visible at the hovered point.
    vil_image_view<vxl_byte> vxl_frame = 
        vil_load(filenames_[vframe_index].c_str());

    // Normalize the gain and offset of greylevels
    const double gain = aligner_.gain(vframe_index);
    const double offset = aligner_.offset(vframe_index);
    if ( (vcl_abs(gain - 1.0) > 1e-6) ||
         (vcl_abs(offset - 0.0) > 1e-6) )
    {
      vil_math_scale_and_offset_values(vxl_frame, gain, offset);
    }

    // Find the mouse position in the coordinate frame of the currently
    // displayed frame.
    const QPoint imagePosition = mousePosition_.toPoint() -
                                 stitcher_->position(vframe_index);

    // Crop the region around the current mouse position.
    vxl_frame = vil_crop(vxl_frame, imagePosition.x()-64, 128, 
                                    imagePosition.y()-64, 128);

    // Update the liveView image.
    ui_.liveView->setPixmap(pixmap_from(vxl_frame, 2.0));
  }

  // Update the frame in the selected view (i.e. showing the region directly
  // around the position where the mouse was last clicked).
  static unsigned sframe_index = 0;

  if (selected_frames_.empty())
    ui_.selectedView->setPixmap(QPixmap());
  else
  {
    // Go to the next visible frame.
    if ( (selected_frames_.front() <= sframe_index) &&
                                     (sframe_index < selected_frames_.back()) )
      ++sframe_index;
    else
      sframe_index = selected_frames_.front();

    // Load the frame visible at the selected point.
    vil_image_view<vxl_byte> vxl_frame = 
        vil_load(filenames_[sframe_index].c_str());

    // Normalize the gain and offset of greylevels
    const double gain = aligner_.gain(sframe_index);
    const double offset = aligner_.offset(sframe_index);
    if ( (vcl_abs(gain - 1.0) > 1e-6) ||
         (vcl_abs(offset - 0.0) > 1e-6) )
    {
      vil_math_scale_and_offset_values(vxl_frame, gain, offset);
    }

    // Find the selected position in the coordinate frame of the currently
    // displayed frame.
    const QPoint imagePosition = selectedPosition_.toPoint() -
                                 stitcher_->position(sframe_index);

    // Crop the region around the currently selected position.
    vxl_frame = vil_crop(vxl_frame, imagePosition.x()-64, 128, 
                                    imagePosition.y()-64, 128);

    // Update the selectedView image.
    ui_.selectedView->setPixmap(pixmap_from(vxl_frame, 2.0));
  }
}
 
//: React to the mouse being moved over the mosaic by changing the point
//  at which images are cropped for the liveView.
void ncm_qstitch_main::on_mosaicView_mouseMoved(QPointF scenePos)
{
  mousePosition_ = scenePos;

  if (stitcher_ != NULL)
  {
    // Get the list of consecutive frames in which a 128x128 box around the
    // current mouse position are all visible.
    visible_frames_ = find_consecutive(stitcher_->images_at(scenePos));

    // Update the status bar.
    if (!visible_frames_.empty())
    {
      ui_.liveLabel->setText("Live View (" + 
                             QString::number(visible_frames_.size()) + 
                             " frames)");
    }
    else
    {
      ui_.liveLabel->setText("Live View: No Frames");
    }
  }
  else
  {
    visible_frames_.clear();
    ui_.liveLabel->setText("Live View: No Frames");
  }
}
 
//: React to the mouse being clicked on the mosaic by changing the point
//  at which images are cropped for the selectedView.
void ncm_qstitch_main::on_mosaicView_mousePressed(QPointF scenePos)
{
  if (stitcher_ != NULL)
  {
    // Delete any existing item and create a new one in the desired place.
    // Because of the weird coordinate frames used, this is easier than moving 
    // the existing item.
    if (flowRectItem_ != NULL)
    {
      scene_.removeItem(flowRectItem_);
      delete flowRectItem_;
    }
    QRectF flowRect(scenePos.x()-64, scenePos.y()-64, 128, 128);
    flowRectItem_ = scene_.addRect(flowRect, QPen(QColor(255,0,0)));

    // Get the list of consecutive frames in which a 128x128 box around the
    // clicked mouse position are all visible.
    selected_frames_ = find_consecutive(stitcher_->images_at(scenePos));
    if (!selected_frames_.empty())
    {
      ui_.selectedLabel->setText("Selected View (" + 
                                 QString::number(selected_frames_.size()) + 
                                 " frames)");

      // Allow flow estimation if the number of visible frames is above 40.
      // (Of course, this is arbitrary and shouldn't be hard-coded.)
      if (selected_frames_.size() > 40)
      {
        ui_.estimateFlow->setEnabled(true);
        selectedPosition_ = scenePos;
      }
      else
      {
        ui_.estimateFlow->setEnabled(false);
        selectedPosition_ = QPointF();
      }
    }
    else
    {
      ui_.selectedLabel->setText("Selected View: No Frames");
      ui_.estimateFlow->setEnabled(false);
      selectedPosition_ = QPointF();
      selected_frames_.clear();
    }
  }
}
 
//: Estimate the optical flow in the region surrounding the selected point.
void ncm_qstitch_main::on_estimateFlow_clicked()
{
  // Shouldn't need to check that selected_frames_ is nonempty.
  // If it isn't, this button should have been disabled.
  // (But I will, just in case.)

  if (selected_frames_.empty())
  {
    ui_.estimateFlow->setEnabled(false);
    return;
  }

  // Create a list of images of the region surrounding the selected point.
  vcl_vector< vil_image_view<vxl_byte> > image_stack;
  for (unsigned i = 0; i < selected_frames_.size(); ++i)
  {
    // Load the original frame from disk.
    unsigned frame = selected_frames_[i];
    vil_image_view<vxl_byte> vxl_frame = vil_load(filenames_[frame].c_str());

    // Normalize the gain and offset of greylevels
    const double gain = aligner_.gain(frame);
    const double offset = aligner_.offset(frame);
    if ( (vcl_abs(gain - 1.0) > 1e-6) ||
         (vcl_abs(offset - 0.0) > 1e-6) )
    {
      vil_math_scale_and_offset_values(vxl_frame, gain, offset);
    }

    // Get the position of the selected point in the image coordinate frame.
    const QPoint imagePosition = selectedPosition_.toPoint() -
                                 stitcher_->position(frame);

    // Add the cropped region to the image_stack.
    vil_image_view<vxl_byte> cropped_frame;
    cropped_frame = vil_crop(vxl_frame, imagePosition.x()-64, 128, 
                                        imagePosition.y()-64, 128);
    image_stack.push_back(cropped_frame);
  }

  // Create a flow field object and estimate the field from the image stack
  // (TODO: Move this to a separate thread so that it doesn't hang the GUI
  //        while it's processing.)
  ncm_flow_field flow_field;
  ncm_flow_estimator_cholmod flow_estimator;
  flow_estimator.estimate_flow_from(image_stack, flow_field, 2);

  // Show the result in the flowView.
  vil_image_view<vxl_byte> flow_image_vxl = flow_field.as_colormap();
  ui_.flowView->setPixmap(pixmap_from(flow_image_vxl, 2.0));
}
 
//: Pause execution of the program so that you can look at the internal
//  variables of the classes.
void ncm_qstitch_main::on_debug_clicked()
{
  int dummy_variable = 0;
}

//
// Private methods
//
 
//: Create the status bar's labels and progress bar.
void ncm_qstitch_main::initializeStatusbar()
{
  statusBar_main_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_main_.setLineWidth(1);
  statusBar_main_.setText("");
  statusBar_main_.setContentsMargins(4, 0, 4, 0);
  statusBar()->addWidget(&statusBar_main_, 1);

  statusBar_progress_.setFixedWidth(160);
  statusBar()->addWidget(&statusBar_progress_);
}
 
//: Update status bar
void ncm_qstitch_main::updateStatus(const QString& status, 
                                    double progress /* in [0..1] */)
{
  statusBar_main_.setText(status);

  // Reset progress bar if negative value supplied.
  if (progress < 0.0)
    statusBar_progress_.reset();
  else
  {
    const int progress_min = statusBar_progress_.minimum();
    const int progress_range = statusBar_progress_.maximum() - progress_min;
    statusBar_progress_.setValue(progress_min + progress*progress_range);
  }

  statusBar()->update();
}
 
//: Set up the thread that contains the aligner object.
void ncm_qstitch_main::initializeAlignerThread()
{
  aligner_thread_.setObjectName("aligner_thread");

  aligner_.moveToThread(&aligner_thread_);

  // Register data types before they are used in signals
  qRegisterMetaType< vcl_vector<vcl_string> >("vcl_vector<vcl_string>");

  // This -> aligner
  QObject::connect( this, SIGNAL(readyToAlign()),
                    &aligner_, SLOT(alignNextFrame()) ); 

  // Aligner -> this
  QObject::connect( &aligner_, 
                      SIGNAL(framesAligned(int, int, double, double, double)),
                    this, 
                      SLOT(onFramesAligned(int, int, double, double, double)) );
  QObject::connect( &aligner_, 
                      SIGNAL(alignmentFinished()),
                    this, 
                      SLOT(onAlignmentFinished()) );
}
 
//: Set up the thread for the stitcher.
void ncm_qstitch_main::initializeStitcherThread()
{
  // Stitcher must be created here as it must be given an image size which
  // isn't known at compile time.
  stitcher_ = new ncm_qmosaic_maker(aligner_.ni(), aligner_.nj());

  stitcher_thread_.setObjectName("mosaic_thread");

  stitcher_->moveToThread(&stitcher_thread_);

  // Register data types so that they can be used in signals.
  qRegisterMetaType< vcl_vector<vcl_string> >("vcl_vector<vcl_string>");
  qRegisterMetaType< vcl_vector<int> >("vcl_vector<QPoint>");

  // this -> stitcher_
  QObject::connect( this, 
                      SIGNAL(makeMosaic(vcl_vector<vcl_string>,
                                        vcl_vector<QPoint>)),
                    stitcher_, 
                      SLOT(makeMosaicFrom(vcl_vector<vcl_string>,
                                          vcl_vector<QPoint>)) ); 

  // stitcher_ -> this
  QObject::connect( stitcher_, SIGNAL(mosaicUpdated()),
                    this, SLOT(onMosaicUpdated()) );
  QObject::connect( stitcher_, SIGNAL(mosaicFinished()),
                    this, SLOT(onMosaicFinished()) );
}
 
//: Add a picture to the scene at a specified (x,y) location.
QGraphicsPixmapItem* ncm_qstitch_main::add_pixmap(
  QPixmap pixmap, 
  QPointF position)
{
  QGraphicsPixmapItem* item = scene_.addPixmap(pixmap);
  item->setPos(position);

  return item;
}
 
//: Add the mosaic image - nothing else - and resize the view such that the
//  image fills the view.
void ncm_qstitch_main::show_as_mosaic()
{
  scene_.clear();

  QGraphicsPixmapItem* p = scene_.addPixmap(pixmap_from(stitcher_->mosaic()));
  scene_.setSceneRect(p->boundingRect());

  ui_.mosaicView->fitInView(scene_.sceneRect(), Qt::KeepAspectRatio	);
}
 
//: Clean up, ready for another batch of images.
void ncm_qstitch_main::reset()
{
  animateTimer_.stop();

  aligner_thread_.exit();
  stitcher_thread_.exit();

  scene_.clear();
  frame_items_.clear();
  border_item_ = NULL;
  flowRectItem_ = NULL;
  pixelCountItem_ = NULL;
  maskItem_ = NULL;

  filenames_.clear();
  visible_frames_.clear();
  selected_frames_.clear();

  mousePosition_ = QPointF();
  selectedPosition_ = QPointF();

  ui_.liveView->setPixmap(QPixmap());
  ui_.selectedView->setPixmap(QPixmap());
  ui_.flowView->setPixmap(QPixmap());

  n_aligned_ = 0;
  n_stitched_ = 0;

  delete stitcher_;
  stitcher_ = NULL;
}
 
//: Find all of the PNG images in a given folder.
void ncm_qstitch_main::get_filenames_from(vcl_string image_file)
{
	// load the markup
	sequence_.t_read(image_file);

	vcl_string image_dir = vul_file::dirname(image_file);

	int num_frames = sequence_.num_frames();
	for (int i = 0; i < num_frames; ++i)
	{
		vcl_string frame_name = image_dir + "/" + sequence_.frame_header(i).frame_name_;
		//vcl_cout << frame_name << vcl_endl;
		filenames_.push_back(frame_name);
	}

  /*vcl_string search_str = image_dir + "/frame*.png";

  filenames_.clear();
  vul_file_iterator fn(search_str.c_str());
  for (fn; fn; ++fn)
    filenames_.push_back(fn());
	*/

  // Ensure filenames are sorted (vul_file_iterator doesn't)
  vcl_sort(filenames_.begin(), filenames_.end());
}
 
//: Enable controls that should only be enabled once a valid folder is 
//  specified
void ncm_qstitch_main::update_controls()
{
  bool enabled;
  if (n_images() > 0)
  {
    enabled = true;
    ui_.frameSlider->setMaximum(n_images()-1);
  }
  else
  {
    enabled = false;
  }

  ui_.frameSlider->setEnabled(enabled);
  ui_.coarseBox->setEnabled(enabled);
}
 
//: Load the motor readings from the associated log file.
bool ncm_qstitch_main::load_motor_data()
{
	/*
  // Get the filename of the sequence_properties file
  vcl_string search_str = root_dir_ + "/sequence_properties*.txt";
  vcl_string sp_filename;
  vul_file_iterator fn(search_str.c_str());
  for (fn; fn; ++fn)
    sp_filename = fn();

  if ( !sp_filename.empty() )
  {
    vcl_ifstream sequence_properties(sp_filename.c_str());

    char s_tmp[256];

    // Discard the header (15 lines)
    for (unsigned i = 0; i < 11; ++i)
      sequence_properties.getline(s_tmp, 256);
		*/

	// Now read the numbers
  timestamps_.clear();
  motor_positions_.clear();
  motor_zoom_.clear();
  sharpness_vec_.clear();

	int num_frames = sequence_.num_frames();
	timestamps_.resize(num_frames);
  motor_positions_.resize(num_frames);
  motor_zoom_.resize(num_frames);
  sharpness_vec_.resize(num_frames);

	// Conversion from mm (the units of the motor positions) and pixels in the
  // scene. This may not be totally accurate.
  const double pixels_per_mm = 1166.1;//803.85;

	// These are needed to discard the first few lines if they are not
  // valid numbers.
  bool is_first_valid_position = true;
  QPointF datum(0,0);

	for (int i = 0; i < num_frames; ++i)
	{
		ncm_video_frame_header &header = sequence_.frame_header(i);

		QPointF pos(header.motor_position_[0],
                header.motor_position_[1]);
			
		bool const position_is_valid = true; //TO DO - what are valid values?

		if (!position_is_valid)
        motor_positions_[i] = QPointF(0.0,0.0);
    else
    {
      // Transform to reflect the reversed motors and differences in scale.
      // (TODO: Apply the negation for reversed motors in the capture 
      //        software, as we can't assume it here and it isn't recorded.)
      pos = -pos * pixels_per_mm;
			pos.setY(-pos.y());

      // Set the datum position to the first valid position.
      if (is_first_valid_position)
      {
        datum = pos;
        is_first_valid_position = false;
      }

      // Store subsequent position relative to datum.
      // Note that the first position will therefore be (0,0).
      motor_positions_[i] = pos - datum;
    }

    // Read the motor positions (zoom).
		motor_zoom_[i] = header.motor_position_[2];

    // Read the estimated sharpness (not used at the moment).
		sharpness_vec_[i] = header.sharpness_;

		//Read the timestamp
		timestamps_[i] = header.frame_time_.toString("hhmmsszzz").toInt();
  }
	return true;
		/*
    while ( !sequence_properties.eof() )
    {
      // Read the first number on the line as a string.
      vcl_string num_str;
      sequence_properties >> num_str;

			vcl_cout << "Frame string: " << num_str << vcl_endl;

      // Convert it to an unsigned (the sample number, which starts at 1).
      unsigned num = vul_string_atoi(num_str);
      
      bool const valid_number = (num != 0);
      if (valid_number)
      {
        // Read the timestamp.
        vcl_string timestamp_str;
        sequence_properties >> timestamp_str;
        timestamps_.push_back(vul_string_atoi(timestamp_str));

        // Read the motor positions (X and Y).
        vcl_string motorX_str;
        sequence_properties >> motorX_str;
        vcl_string motorY_str;
        sequence_properties >> motorY_str;
        QPointF pos(vul_string_atof(motorX_str),
                    vul_string_atof(motorY_str));

        bool const position_is_valid = ((motorX_str != "-1.#IND") && 
                                        (motorY_str != "-1.#IND"));

        if (!position_is_valid)
          motor_positions_.push_back( QPointF(0.0,0.0) );
        else
        {
          // Transform to reflect the reversed motors and differences in scale.
          // (TODO: Apply the negation for reversed motors in the capture 
          //        software, as we can't assume it here and it isn't recorded.)
          pos = -pos * pixels_per_mm;

          // Set the datum position to the first valid position.
          if (is_first_valid_position)
          {
            datum = pos;
            is_first_valid_position = false;
          }

          // Store subsequent position relative to datum.
          // Note that the first position will therefore be (0,0).
          motor_positions_.push_back(pos - datum);
        }

        // Read the motor positions (zoom).
        vcl_string motorZ_str;
        sequence_properties >> motorZ_str;
        motor_zoom_.push_back(vul_string_atof(motorZ_str));

        // Read the estimated sharpness (not used at the moment).
        vcl_string sharpness_str;
        sequence_properties >> sharpness_str;
        sharpness_vec_.push_back(vul_string_atof(sharpness_str));
      }
    }
    sequence_properties.close();

    // Return true if we read all the values successfully.
    return true;
  }

  // Return false if there was a problem.
  return false;*/
}
 
//: Return the number of images in the sequence. Currently, this is read from
//  the UI spinbox but could be redefined.
unsigned ncm_qstitch_main::n_images()
{
  return ui_.nFramesSpin->value();
}
 
//: Draw a path in the view of the read motor positions.
void ncm_qstitch_main::draw_motor_positions()
{
  if (motor_positions_.size() > 0)
  {
    QPainterPath path;

    path.moveTo(motor_positions_[0]);
    for (unsigned i = 1; i < motor_positions_.size(); ++i)
      path.lineTo(motor_positions_[i]);

    scene_.addPath( path, QPen(QColor(255,0,0)) );
  }
}
 
//: Display a directory selection box with which the user will select the input
//  image folder.
void ncm_qstitch_main::select_folder()
{

	QString selfilter;
	QString folderName = QFileDialog::getOpenFileName(
    this, "Load sequence", "C:/isbe/nailfold/camera_capture",
    tr("Sequence files (*.txt)" ), // filter
    &selfilter );

  /*// get filename from dialog box
	QString folderName = QFileDialog::getExistingDirectory(
                          this, 
                          tr("Select folder to save images in"), 
                          "C:/isbe/nailfold/camera_capture", 
                          QFileDialog::ShowDirsOnly);*/


  // Do nothing more if the box was cancelled.
  if (folderName.isNull())
    return;
	
  reset();

  root_dir_ = folderName.toStdString();
  get_filenames_from(root_dir_);

  if (filenames_.size() == 0)
  {
    // Do something if there are no images in the folder.
    // (TODO: Probably should be more done here.)
    ui_.nFramesSpin->setEnabled(false);
    vcl_cout << "No images to process" << vcl_endl;
    return;
  }
  else
  {
    // Initialize the aligner with filenames
    aligner_.setFilenames(filenames_);

    // Set up some of the GUI controls to match the sequence length.
    ui_.nFramesSpin->setMaximum(filenames_.size());
    ui_.nFramesSpin->setValue(filenames_.size());
    ui_.nFramesSpin->setEnabled(true);

    update_controls();
    load_motor_data();
    draw_motor_positions();

    mosaicRect_ = scene_.sceneRect();
    ui_.mosaicView->fitInView(scene_.sceneRect(), Qt::KeepAspectRatio);
  }
}
 
//: Return a pixmap object from the input image, scaled in size by scale, and 
//  with grey levels adjusted by gain and offset.
QPixmap ncm_qstitch_main::pixmap_from(
  vil_image_view<vxl_byte> image_in,
  double scale /* = 1.0 */,
  double gain /* = 1.0 */,
  double offset /* = 0.0 */)
{
  vil_image_view<vxl_byte> vxl_image;

  if ( (vcl_abs(gain - 1.0) > 1e-6) ||
       (vcl_abs(offset - 0.0) > 1e-6) )
  {
    // Scale and offset a copy of the input image, since a vil_image_view acts
    // as a pointer and will change the aligner's copy of the image.
    vxl_image.deep_copy(image_in);
    vil_math_scale_and_offset_values(vxl_image, gain, offset);
  }
  else
    vxl_image = image_in;

  QImage qt_image;
  qcore_convert_image(qt_image, vxl_image);
  
  QPixmap pixmap;
  int scaled_width = scale * qt_image.width();
  if (scaled_width == qt_image.width())
    pixmap.convertFromImage(qt_image);
  else
    pixmap.convertFromImage(qt_image.scaledToWidth(scaled_width));

  return pixmap;
}
 
//: Move the frame items to their correct positions.
void ncm_qstitch_main::update_frameitem_positions(
  unsigned first /* = 0 */)
{
  const unsigned last = frame_items_.size();

  QPoint position(0,0);
  for (unsigned i = 0; i < last; ++i)
  {
    if ( (i >= first) && (frame_items_[i] != NULL) )
      frame_items_[i]->setPos(position);

    position += aligner_.displacement(i);
  }
}
 
//: Begin the alignment process on the sequence.
void ncm_qstitch_main::start_alignment()
{
  aligner_.reset();

  scene_.clear();
  //scene_.setSceneRect(0,0, 1,1);
  scene_.setSceneRect(0,0, 0,0);
  
	frame_items_.clear();
  frame_items_.resize(n_images(), NULL);
  QPixmap pixmap = pixmap_from(aligner_.destination_image());
  frame_items_[0] = add_pixmap(pixmap, QPoint(0,0));

  // Set the main view to include a very large rectangle.
  ui_.mosaicView->fitInView(QRectF(-2000,1000,4000,2000), Qt::KeepAspectRatio);

  // Start the clock and align the first frames. Subsequent pairs of frames 
  // will be aligned when the signal is received from the aligner.
  start_time = vcl_clock();
  aligner_thread_.start();
  emit readyToAlign();
}