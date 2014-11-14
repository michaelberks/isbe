#ifndef NCM_QCAPTURE2_GUI_H
#define NCM_QCAPTURE2_GUI_H

#include <QMainWindow>
#include <QWheelEvent>
#include <QTime>
#include <QFileDialog>
#include <QBitmap>
#include <QGraphicsPixmapItem>
#include <QMessageBox>
#include <QDebug>
#include <QLabel>
#include <QProgressBar>

#include <vcl_iostream.h>
#include <vcl_cmath.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vcl_vector.h>
#include <vcl_cstddef.h>

#include <vil/vil_image_view.h>

#include <vsl/vsl_indent.h>
#include <vsl/vsl_binary_io.h>

#include <vul/vul_string.h>

#include <nailfold/ncm_registry_accessor.h>

#include <nailfold/qtools/ncm_qapt_server.h>
#include <nailfold/qtools/ncm_qframe_aligner.h>

#include <nailfold/qtools/ncm_cameras/dmk_listener.h>
#include <nailfold/qtools/ncm_cameras/ncm_dmk_camera.h>

#include "ncm_qcapture2_saver.h"
#include "ncm_qcapture2_processor.h"
#include "ncm_qcapture2_data_manager.h"

#include "ui_ncm_qcapture2_gui.h"


class ncm_qcapture2_gui : public QMainWindow
{
  Q_OBJECT

public: // Methods
  
  ncm_qcapture2_gui(QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture2_gui();

  //: Non UI members manipulated by the UI
  void get_frame();

protected:
  
  //void timerEvent( QTimerEvent * );
  void closeEvent(QCloseEvent* ev);

signals:

  void frame_tagged(int, bool, bool, bool);
  void frame_logged(int);

  void apt_move_zero(ncm_qapt_server::apt_axis, bool);
  void apt_move_home(ncm_qapt_server::apt_axis, float, bool);
  void apt_move_home(ncm_qapt_server::apt_axis, bool);
  void apt_move_to(ncm_qapt_server::apt_axis, float, bool);
  void apt_move_by(ncm_qapt_server::apt_axis, float, bool);
  
private slots:

	//Load a camera
	void on_loadCameraButton1_clicked();
	void on_loadCameraButton2_clicked();

	//Change which camera is active
	void change_active_camera( int cam_num );

	//Display controls
	void on_brightnessSlider_sliderMoved(int value);
	void on_contrastSlider_sliderMoved(int value);
	void on_stopLiveButton_clicked();
	void on_fliphImage1Button_toggled();
	void on_fliphImage2Button_toggled();
	void on_flipvImage1Button_toggled();
	void on_flipvImage2Button_toggled();

	//Camera controls
	void on_frameRateSelection_currentIndexChanged(const QString &text);
	void on_gainSpinBox_valueChanged(double value);
	void on_gainSlider_sliderMoved(int value);
	void on_autoGain_stateChanged(int state);
	void on_exposureSpinBox_valueChanged(double value);
	void on_exposureSlider_sliderMoved(int value);
	void on_autoExposure_stateChanged(int state);

	//Save sequence controls
	void on_sequenceLengthSpinBox_valueChanged(double value);
	void on_saveButton_clicked();
	void on_saveDirSelect_clicked();
	void on_sequenceNameTextEdit_textChanged();
	void on_patientNameTextEdit_textChanged();
	void on_handComboBox_currentIndexChanged (const QString &text);
	void on_digitComboBox_currentIndexChanged (const QString &text);
	void on_magComboBox_currentIndexChanged (const QString &text);
	void on_camera1SuffixTextEdit_textChanged();
	void on_camera2SuffixTextEdit_textChanged();
	void on_useCamera1Button_clicked();
	void on_useCamera2Button_clicked();

	//Process new frame grabbed from the camera
	void redraw_scene( int cam_num );
	void redraw_scene();
	void draw_diff_frame();
	void on_processing_finished();

  void add_properties_header();
  void log_frame_properties(int frame_number);
  void tag_frame(int cam_num, bool fliph, bool flipv);
	void tag_frame1();
	void tag_frame2();
  void on_frame_saved(int frame_number, int cam_num);
  void on_saving_finished(int cam_num);

  //: Slots for motor controls
  void onAptTimer();
  void on_aptXCombobox_activated(int index);
  void on_aptYCombobox_activated(int index);
  void on_aptZCombobox_activated(int index);

  void on_aptXReverse_stateChanged(int state);
  void on_aptYReverse_stateChanged(int state);
  void on_aptZReverse_stateChanged(int state);

  void on_aptXVelocity_sliderReleased();
  void on_aptYVelocity_sliderReleased();
  void on_aptZVelocity_sliderReleased();

  void on_aptXHome_clicked();
  void on_aptYHome_clicked();
  void on_aptZHome_clicked();

  // Connect to motors.
  void on_connectToMotors_clicked();

  // Initiate the calibration motion.
  void on_calibrate_clicked();

  // Autofocus on scene.
  void on_autoFocus_clicked();
  void onAptMoveComplete(ncm_qapt_server::apt_axis);
  void onLowerSharpnessThresholdCrossed();
  void onUpperSharpnessThresholdCrossed();

	//Slots for the difference image
	void on_computeDiffButton_clicked();
	void on_diffBrightnessSlider_sliderMoved(int value);
	void on_diffContrastSlider_sliderMoved(int value);

  //: Other slots
  void on_graphicsView_clicked(double x, double y);
  void on_graphicsView_zoomed(int delta);

private: // Methods

	//load a camera
	bool load_camera( int cam_num );

  //: Current state of the autofocus procedure
  enum autoFocus_state_enum {
    autoFocus_InvalidFirst = -1,
    autoFocus_Inactive,
    autoFocus_GetDirectionForward,
    autoFocus_GetDirectionReverse,
    autoFocus_FindPeak,
    autoFocus_GoToPeak,
    autoFocus_Error,

    autoFocus_First = autoFocus_Inactive,
    autoFocus_Last = autoFocus_Error
  };

  //: Current state of the calibration procedure.
  enum calibrate_state_enum {
    calibrate_InvalidFirst = -1,
    calibrate_Inactive,
    calibrate_Up,
    calibrate_Right,
    calibrate_Down,
    calibrate_Left,
    calibrate_Error,

    calibrate_First = calibrate_Inactive,
    calibrate_Last = calibrate_Error
  };

  void update_framerate_controls();
  void update_exposure_controls();
  void update_gain_controls();

  //: Create status bar components.
  void initializeStatusbar();

  //: Update text and progress.
  void updateStatus(
    const QString& status, 
    double progress = -1);

  void updateProgress(
    double progress);

  //: Initialize aligner thread.
  void initializeAlignerThread();

  //: Update raw colormap to enhance contrast
	void update_contrasts(int brightness, int contrast,  int & c_min, int & c_max);
  void update_raw_colour_table(int brightness, int contrast, QVector<QRgb> &raw_colour_table);
	void update_diff_colour_table(int brightness, int contrast);

	//Write out text file of camera properties
	void write_properties( vcl_ostream& tfs );
	void write_camera_properties( vcl_ostream& tfs , int cam_num );

  void draw_pixmap(QImage* qimage, int cam_num);
  void draw_crosshairs(int cam_num);

  void initialize_registry_accessor();
  void initialize_saver();
  void initialize_apt_timer();
  void initialize_apt_spinboxes();
  void update_apt_comboboxes();

  void set_apt_targetX(double x);
  void set_apt_targetY(double y);
  void set_apt_targetZ(double z);

  void autoFocus(autoFocus_state_enum state);
  void handleAutoFocusState();
  void handleCalibrateState();

  void initialize_processor_thread();
  void initialize_saver_thread();

  void get_apt_registry_values();
  void initialize_apt_thread();

  void connect_signals_to_slots();

private: // Variables

  //: User interface object
  Ui::NailfoldMainWindow ui;

  //: Status bar components.
  QLabel statusBar_main_;
  QProgressBar statusBar_progress_;

	//: Cursor graphic to use when marking up
  QCursor markup_cursor_;

	//: pixmap item for raw image in scene
  QVector<QGraphicsPixmapItem*> raw_pixmap_items_;

	//: pixmap item for raw image in scene
  QGraphicsPixmapItem* raw_pixmap_item_d_;

  //: Colormap for raw images
  QVector < QVector<QRgb> > raw_colour_tables_;
	QVector<QRgb> raw_colour_table_d_;

  //: canvas for each camera on which to draw stuff
  QVector<QGraphicsScene*> scenes_;

	//: canvas to draw difference between 2 camera views
	QGraphicsScene scene_d_;

	//: brightnessa nd contrast values for each camera
	QVector<int> brightness_;
	QVector<int> contrast_;

	//Camera objects
	QVector<ncm_dmk_camera*> cameras_;

	//: Index to current camera
	int current_camera_;

	//: Flags showing which cameras to save from
	QVector<bool> save_from_camera_;

	//: Flags showing which cameras have diff frames tagged
	QVector<bool> diff_frame_;

	//: Flags showing which cameras have been flipped
	QVector<bool> flipped_h_;
	QVector<bool> flipped_v_;

	//Saver objects
	QVector<ncm_qcapture2_saver*> savers_;

	//Saver thread
	QThread saver_thread_;

	//Processor object
	ncm_qcapture2_processor processor_;

	//Processor thread
	QThread processor_thread_;

  //: Constants related to motor control
  const bool motors_are_live_;
  const double pixels_per_mm_;

  //: Motor controller
  //ncm_qapt_server apt_;

	//Processor thread
	QThread apt_thread_;

  //: Timer for updating motor position feedback
  QTimer apt_timer_;

  //: Log file for focus values
  vcl_ofstream focus_log_;

  //: Index of Finite State Machine that controls autofocus
  autoFocus_state_enum autoFocus_state_;

  //: Index of Finite State Machine that controls autofocus
  calibrate_state_enum calibrate_state_;

  //: Position of point where sharpness was maximized.
  float autoFocus_peak_position_;

  //: Number of samples to use to compute running average of sharpness
  const unsigned n_sharpness_;

  //: Registry accessor
  ncm_registry_accessor ra_;

  //: Image of dirt on the lens
  vil_image_view<int> artefact_image_;

  //: Aligner for registering frames
  ncm_qframe_aligner aligner_;

	//: Aligner thread
	QThread aligner_thread_;

  //: File handle for properties
  vcl_ofstream property_fs_;

	//Widgets added manually
	QButtonGroup* currentCameraRadioGroup;
};

#endif // NCM_QCAPTURE2_GUI_H
