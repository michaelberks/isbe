#ifndef NCM_QSTITCH_MAIN_H
#define NCM_QSTITCH_MAIN_H

#include "ui_ncm_qstitch_main.h"

#include <vcl_ctime.h>

#include <vbl/vbl_shared_pointer.h>

#include <vgl/vgl_point_2d.h>

#include <nailfold/qtools/ncm_qframe_aligner_file.h>
#include <nailfold/qtools/ncm_qmosaic_maker.h>
#include <nailfold/qtools/ncm_image_sequence.h>

#include <QMainWindow>
#include <QLabel>
#include <QProgressBar>
#include <QThread>
#include <QGraphicsScene>
#include <QTimer>

class ncm_qstitch_main : public QMainWindow
{
  Q_OBJECT

public: // Methods

  //: Default constructor
  ncm_qstitch_main(
    QWidget *parent = 0, Qt::WFlags flags = 0);

  //: Destructor
  ~ncm_qstitch_main();

signals:
  
  //: Signal emitted when the next frame should be aligned with the sequence.
  void readyToAlign();

  //: Signal emitted when the stitching should start.
  void makeMosaic(
    vcl_vector<vcl_string>,
    vcl_vector<QPoint>);

public slots:

protected: // Methods

private slots:

  //: Select the input image folder.
  void on_actionSelectFolder_triggered();

  //: Select the input image folder.
  void on_selectFolder_clicked();

  //: Limit the number of frames aligned and stitched into a mosiac.
  void on_nFramesSpin_valueChanged(int nFrames);

  //: Erase existing displacements in the aligner class.
  void on_clearDisplacements_clicked();

  //: Create a mask so that dirt on the lens (which is in the same place at every
  //  frame and can therefore be detected by averaging all the frames) is masked
  //  out before alignment.
  void on_createMask_clicked();
   
  //: Show/hide the mask in the scene view.
  void on_showMask_toggled(bool checked);
  
  //: Copy the displacements from the motor positions to the aligner as an
  //  initialization. In practice, 'good' sequences don't need this, though
  //  sequences with large regions of uniform texture may benefit from it.
  void on_useAllMotorPositions_clicked();
  
  //: Change the maximum distance allowed between consecutive frames.
  void on_radiusSpin_valueChanged(int value);

  //: Align the frames in the sequence using the fast Sobel method on heavily
  //  downsampled images.
  void on_coarseAlign_clicked();

  //: Align the frames in the sequence using the slower Gaussian 1st derivatives
  //  method on non-downsampled images.
  void on_fineAlign_clicked();

  //: React to two frames having been aligned by adding the more recent image
  //  to the scene at its aligned position.
  void onFramesAligned(
    int di, int dj, 
    double scale, double offset, 
    double mse);

  //: React to all frames having been aligned.
  void onAlignmentFinished();

  //: Change the horizontal displacement of the currently selected frame (which
  //  also affects all subsequent frames).
  void on_horizDispSpin_valueChanged(int dx);
 
  //: Change the vertical displacement of the currently selected frame (which
  //  also affects all subsequent frames).
  void on_vertDispSpin_valueChanged(int dy);

  //: Copy the displacement from the motor positions for this frame only.
  void on_useMotorPositions_clicked();

  //: Change the currently selected frame in response to sliding the slider.
  void on_frameSlider_valueChanged(int frame);

  //: Create a mosaic of the aligned images by blending.
  void on_blendButton_clicked();

  //: React to the mosaic being updated.
  void onMosaicUpdated();

  //: React to the mosaic being finished by displaying in the scene view.
  void onMosaicFinished();

  //: Show the pixel count image (i.e. a monochrome in which pixels visible in
  //  the most number of frames are shown in white, whereas those visible in the
  //  least are shown in black).
  void on_showPixelCount_toggled(bool checked);

  //: Save the mosaic to a desired filename.
  void on_saveImage_clicked();

  //: Estimate the optical flow in the region surrounding the selected point.
  void on_estimateFlow_clicked();
  
  //: React to the timer timeout by updating the frame in the respective views.
  void onAnimateTimer();

  //: React to the mouse being moved over the mosaic by changing the point
  //  at which images are cropped for the liveView.
  void on_mosaicView_mouseMoved(QPointF scenePos);

  //: React to the mouse being clicked on the mosaic by changing the point
  //  at which images are cropped for the selectedView.
  void on_mosaicView_mousePressed(QPointF scenePos);

  //: Pause execution of the program so that you can look at the internal
  //  variables of the classes.
  void on_debug_clicked();
  
private: // Methods

  //: Create the status bar's labels and progress bar.
  void initializeStatusbar();

  //: Set up the thread that contains the aligner object.
  void initializeAlignerThread();

  //: Set up the thread for the stitcher.
  void initializeStitcherThread();

  //: Update status bar
  void updateStatus(
    const QString& status, 
    double progress /* in [0..1] */);

  //: Add a picture to the scene at a specified (x,y) location.
  QGraphicsPixmapItem* add_pixmap(
    QPixmap pixmap, 
    QPointF position);

  //: Add the mosaic image - nothing else - and resize the view such that the
  //  image fills the view.
  void show_as_mosaic();

  //: Clean up, ready for another batch of images.
  void reset();

  //: Find all of the PNG images in a given folder.
  void get_filenames_from(
    vcl_string image_dir);

  //: Draw a path in the view of the read motor positions.
  void draw_motor_positions();

  //: Enable controls that should only be enabled once a valid folder is 
  //  specified
  void update_controls();

  //: Load the motor readings from the associated log file.
  bool load_motor_data();

  //: Return the number of images in the sequence. Currently, this is read from
  //  the UI spinbox but could be redefined.
  unsigned n_images();

  //: Display a directory selection box with which the user will select the input
  //  image folder.
  void select_folder();

  //: Return a pixmap object from the input image, scaled in size by scale, and 
  //  with grey levels adjusted by gain and offset.
  QPixmap pixmap_from(
    vil_image_view<vxl_byte> vxl_image,
    double scale = 1.0,
    double gain = 1.0,
    double offset = 0.0);

  //: Move the frame items to their correct positions.
  void update_frameitem_positions(
    unsigned first = 0);

  //: Begin the alignment process on the sequence.
  void start_alignment();

private: // Variables

  //: User interface class.
  Ui::MainWindow ui_;

  //: The scene viewed by the viewer.
  QGraphicsScene scene_;

  //: List of all frames displayed in the scene view.
  vcl_vector<QGraphicsItem*> frame_items_;

  //: Red rectangle surrounding the currently selected frame.
  QGraphicsRectItem* border_item_;

  //: Pixmap item displayed in the flow view.
  QGraphicsRectItem* flowRectItem_;

  //: Pixmap item of the pixel count  displayed in the scene view.
  QGraphicsPixmapItem* pixelCountItem_;

  //: Pixmap item of the mask image displayed in the scene view.
  QGraphicsPixmapItem* maskItem_;

  //: Labels to go in status bar. The status bar will contain one main label,
  //  plus two panels at the right hand side to indicate the function of the 
  //  left and right mouse buttons at any given time
  QLabel statusBar_main_;
  QProgressBar statusBar_progress_;

  //: Timer for animation
  QTimer animateTimer_;

  //: Bounding box of the mosaic.
  QRectF mosaicRect_;

  //: The current mouse position.
  QPointF mousePosition_;

  //: The mouse position that was clicked (selected).
  QPointF selectedPosition_;

  //: Root folder for images.
  vcl_string root_dir_;

	//: Image sequence to be aligned
	ncm_image_sequence sequence_;

  //: Filenames of the images to be stitched.
  vcl_vector<vcl_string> filenames_;
  
  //: Coarse displacements of images.
  vcl_vector<QPoint> coarse_displacements_;

  //: Coarse positions of images (cumulative displacements).
  vcl_vector<QPoint> coarse_positions_;

  //: Fine displacements of images.
  vcl_vector<QPoint> fine_displacements_;

  //: Fine positions of images (cumulative displacements).
  vcl_vector<QPoint> fine_positions_;

  //: List of frames visible at the current mouse position.
  vcl_vector<unsigned> visible_frames_;

  //: List of frames visible at the clicked mouse position.
  vcl_vector<unsigned> selected_frames_;

  //: Current position where frames are drawn
  QPoint current_position_;

  //: Aligner QObject
  ncm_qframe_aligner_file aligner_;

  //: Thread for doing the alignment.
  QThread aligner_thread_;

  //: Number of images aligned.
  unsigned n_aligned_;

  //: Mosaic maker QObject
  ncm_qmosaic_maker* stitcher_;

  //: Thread for doing the stitching.
  QThread stitcher_thread_;

  //: Number of images added to the mosaic.
  unsigned n_stitched_;

  //: Timestamps read in from the motor positions log.
  vcl_vector<long> timestamps_;

  //: (x.y) position of the motors read in from the motor positions log.
  vcl_vector<QPointF> motor_positions_;

  //: Z position of the motors, read in from the motor positions log.
  vcl_vector<double> motor_zoom_;

  //: Image sharpness read in from the motor positions log.
  vcl_vector<double> sharpness_vec_;

  //: Start time used when timing operations.
  vcl_clock_t start_time;
};

#endif // NCM_QSTITCH_MAIN_H
