#ifndef NCM_QCAPTURE_GUI_H
#define NCM_QCAPTURE_GUI_H

#include <QMainWindow>
#include <QWheelEvent>
#include <QTime>
#include <QFileDialog>
#include <QBitmap>
#include <QGraphicsPixmapItem>
#include <QMessageBox>
#include <QDebug>
#include <QLabel>
#include <QIcon>
#include <QProgressBar>
#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlQueryModel>
#include <QSortFilterProxyModel>

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

#include <nailfold/qtools/ncm_qapt_server.h>

#include <nailfold/qtools/ncm_cameras/dmk_listener.h>
#include <nailfold/qtools/ncm_cameras/ncm_dmk_camera.h>
#include <nailfold/qtools/ncm_subject.h>
#include <nailfold/qtools/ncm_user.h>
#include <nailfold/ncm_annotation.h>
#include <nailfold/qtools/ncm_qscene.h>
#include <nailfold/qtools/ncm_qimagehandler.h>
#include <nailfold/qtools/ncm_qframe_aligner_file.h>

#include "ncm_qcapture_preferences.h"
#include "ncm_qcapture_saver.h"
#include "ncm_qcapture_processor.h"
#include "ncm_qcapture_aligner.h"
#include "ncm_qcapture_mosaic_maker.h"
#include "ncm_qcapture_data_manager.h"
#include "ncm_qcapture_scene.h"

#include "ui_ncm_qcapture_gui.h"

class ncm_qcapture_gui : public QMainWindow
{
  Q_OBJECT

public: // Methods
  
  ncm_qcapture_gui(ncm_qcapture_preferences&, 
		QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_gui();

  //: Non UI members manipulated by the UI
  void get_frame();

protected:
  
  //void timerEvent( QTimerEvent * );
  void closeEvent(QCloseEvent* ev);
	bool winEvent(MSG * message, long * result);

signals:

  void frame_tagged();
  void frame_logged(int);

  void apt_move_zero(ncm_qapt_server::apt_axis, bool);
  void apt_move_home(ncm_qapt_server::apt_axis, float, bool);
  void apt_move_home(ncm_qapt_server::apt_axis, bool);
  void apt_move_to(ncm_qapt_server::apt_axis, float, bool);
  void apt_move_by(ncm_qapt_server::apt_axis, float, bool);

	void aptXHome_clicked();
	void aptYHome_clicked();
	void aptZHome_clicked();

	void frame_to_align();

	void make_mosaic(vcl_vector<QPoint>);

	void finalize_sequence();

	void review_session_selected(QModelIndex);
	void review_sequence_selected(QModelIndex);
  
private slots:
	void on_actionExit_triggered();

	//Log in/out to the system
	void on_logInPushButton_clicked();
	void on_logOffPushButton_clicked();

	//File menu actions
  void on_actionLoad_camera_triggered();
	void on_actionChange_camera_triggered();
	void on_actionPreferences_triggered();

	// Welcome page control callbacks
	void on_newSessionPushButton_clicked();
	void on_processPushButton_clicked();
	void on_reviewPushButton_clicked();

	//Image session control callbacks
	void on_startRecordingPushButton_clicked();
	void on_stopRecordingPushButton_clicked();
	void on_endSessionPushButton_clicked();

	//Review session calbacks
	void on_finishedReviewButton_clicked();

	//Display controls
	void on_brightnessSlider_sliderMoved(int value);
	void on_contrastSlider_sliderMoved(int value);
	void on_stopLiveButton_clicked();

	// Brightness and contrast controls for reviewing
  void updateContrastControls();
	void on_autoContrast_toggled(bool checked);
  void on_brightnessRSlider_sliderMoved(int value);
  void on_contrastRSlider_sliderMoved(int value);

	//Zoom controls for reviewing
	// Zoom controls
  void updateZoomControls();
  void updateZoomLimits();
  void on_zoomSlider_sliderMoved(int value);
  void on_zoomEdit_editingFinished();
  void on_zoomOutButton_clicked();
  void on_zoomInButton_clicked();
  void applyUserZoom();
  void applyVesselZoom(int vessel_index);

	//Camera controls
	void on_frameRateSelection_currentIndexChanged(const QString &text);
	void on_gainSpinBox_valueChanged(double value);
	void on_gainSlider_sliderMoved(int value);
	void on_autoGain_stateChanged(int state);
	void on_exposureSpinBox_valueChanged(double value);
	void on_exposureSlider_sliderMoved(int value);
	void on_autoExposure_stateChanged(int state);


	//Process new frame grabbed from the camera
	void redraw_scene();
	void draw_mosaic();
	void onProcessingFinished(int n_processed);

  void add_properties_header();
	void log_frame_properties(ncm_video_frame_header header);
  void tag_frame();
  void onFrameSaved(int frame_number);
  void onSavingFinished();
	void onFinalizeSequence();

  void onSceneJoystickMoved(double, double);
  void onSceneJoystickReleased();
  void onSceneDoubleClicked(double, double, Qt::MouseButton);

  //: Slots for motor controls
  void onAptTimer();
  void on_aptXCombobox_activated(int index);
  void on_aptYCombobox_activated(int index);
  void on_aptZCombobox_activated(int index);

  void on_aptXReverse_toggled(bool checked);
  void on_aptYReverse_toggled(bool checked);
  void on_aptZReverse_toggled(bool checked);

  void on_aptXHome_clicked();
  void on_aptYHome_clicked();
  void on_aptZHome_clicked();

  //void on_aptXVelocity_sliderReleased();
	//void on_aptYVelocity_sliderReleased();
	//void on_aptZVelocity_sliderReleased();

  // Connect to motors.
  void on_connectToMotors_clicked();

  // Initiate the calibration motion.
  void on_calibrate_clicked();

  // Autofocus on scene.
  void on_autoFocus_clicked();
  void onAptMoveComplete(ncm_qapt_server::apt_axis);
  void onLowerSharpnessThresholdCrossed();
  void onUpperSharpnessThresholdCrossed();

	//Align frames for mosaic
	void onFrameToAlign(int frame_num);
	void onFramesAligned(
    int di, int dj, 
    double scale, double offset, 
    double mse);
  void onAlignmentFinished();
	void onMosaicUpdated(int frame_num);
	void onMosaicFinished(bool success);
  
	//: Other slots
  //void on_graphicsView_clicked(double x, double y);
  void on_graphicsView_zoomed(int delta);
	//void keyPressEvent(QKeyEvent* event);
	bool eventFilter(QObject *object, QEvent *event);

	//Review slots
	void onSessionSelectionChanged();
	void onSequenceSelectionChanged();
	void on_prevSessionButton_clicked();
	void on_nextSessionButton_clicked();
	void on_prevSequenceButton_clicked();
	void on_nextSequenceButton_clicked();

private: // Methods

	//:Current activity being performed
	enum activity_state_enum {
		activity_Login,
		activity_Home,
		activity_Session,
		activity_Recording,
		activity_RecordingPaused,
		activity_Saving,
		activity_Reviewing
	};
		
  //: Current state of the autofocus procedure
  enum autoFocus_state_enum {
    autoFocus_InvalidFirst = -1,
    autoFocus_Inactive,
    autoFocus_GetDirectionForward,
    autoFocus_GetDirectionReverse,
    autoFocus_FindPeak,
    autoFocus_GoToPeak,
    autoFocus_Interrupt,
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

	void log_in();
	void log_off();
	void start_image_session();
	void end_image_session();

	void start_review_session();
	void start_offline_processing();
	void end_review_session();

	bool load_camera(bool use_existing);
	void disconnect_camera();
	bool connect_to_database(const QString & username, const QString & password);
	void disconnect_database();
	bool load_user(const QString & username);

	void update_session_controls();
	void update_tab_controls();
  void update_framerate_controls();
  void update_exposure_controls();
  void update_gain_controls();
  void disable_motor_controls();

  //: Create status bar components.
  void initializeStatusbar();

  //: Update text and progress.
  void updateStatus(
    const QString& status, 
    double progress = -1);

  void updateProgress(
    double progress);

  //: Initialize aligner thread.
  void initialize_aligner_thread();

  //: Decide contrast limits from slider bar positions.
  void set_contrast_limits(int contrast, int brightness);

  //: Update raw colormap to enhance contrast
  void update_raw_colour_table();

	//Write out text file of camera properties
	void write_properties( vcl_ostream& tfs );

	void initialize_scene();
	void initialize_motors();
  void initialize_saver();
  void initialize_apt_timer();
  void update_apt_comboboxes();

  void set_apt_targetX(double x);
  void set_apt_targetY(double y);
  void set_apt_targetZ(double z);

	void home_x();
  void home_y();
  void home_z();

  void autoFocus(autoFocus_state_enum state);
  void handleAutoFocusState();
  void handleCalibrateState();

  void initialize_processor_thread();
  void initialize_saver_thread();

  void get_apt_sql_values();
  void initialize_apt_thread();

	void update_subject_display();
	void update_session_display();
	void update_review_display();
	void update_review_controls();

	void sync_sequence_name(ncm_image_sequence &sequence);
	void interpolate_motor_positions(ncm_image_sequence &sequence);

  void connect_session_signals_to_slots();
	void disconnect_session_signals_from_slots();

	void connect_review_signals_to_slots();
	void connectSetFocus();
	void disconnect_review_signals_from_slots();

	//: Set the graphicsView zoom factor from a percentage value
  void setZoomFrom(int value);

	void align_full_mosaic(ncm_image_sequence & sequence);
	void make_full_mosaic(ncm_image_sequence &sequence);

private: // Variables

  //: User interface object
  Ui::NailfoldMainWindow ui;

  //: Status bar components.
  QLabel statusBar_main_;
  QProgressBar statusBar_progress_;

	//: Options for life, the universe and everything
	ncm_qcapture_preferences &preferences_;

  float velocityX_;
  float velocityY_;
	float velocityZ_;
    
  //: Canvas on which to draw the camera view
  ncm_qcapture_scene scene_;

	//: Canvas on which to show mosaic
  QGraphicsScene mosaic_scene_;

	//: Canvas on which to show motor positions
  QGraphicsScene motors_scene_;

	//: Green rectangle surrounding the currently selected frame.
  QGraphicsRectItem* current_frame_position_;

	//: Blue rectangles that mark the ath of a mosaic
	vcl_vector<QGraphicsRectItem*> mosaic_frame_positions_;

	//: Position to stick the next mosaic frame during live mosaicing
	QPoint current_mosaic_displacement_; 

	//QGraphicsEllipseItem* z_upper_msg_;
	//QGraphicsEllipseItem* z_lower_msg_;

	//Camera object
	ncm_dmk_camera camera_;

	//Saver object
	ncm_qcapture_saver saver_;

	//Processor object
	ncm_qcapture_processor processor_;
	
	//Saver thread
	QThread saver_thread_;

	//Processor thread
	QThread apt_thread_;

  //Processor thread
	QThread processor_thread_;

	//: Aligner thread
	QThread aligner_thread_;

	//: Offline processing thread
	QThread offline_thread_;

  //: Constants related to motor control
  bool motors_are_live_;

  //: Motor controller
  ncm_qapt_server apt_;

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

	//: Vectors to store motor positions and times during recording (used
	// to correct the time lag between tagging frames with time and position)
	vcl_vector<vcl_vector <double> > motor_positions_;
	vcl_vector<int> motor_times_;
	QTime motor_time0_;

  //: Image of dirt on the lens
  vil_image_view<int> artefact_image_;

  //: Aligner for registering frames
  ncm_qcapture_aligner aligner_;
	ncm_qframe_aligner_file offline_aligner_;

	//: Stitcher for making mosaics (also runs in the aligner thread) 
	ncm_qcapture_mosaic_maker *stitcher_;

	//User of this session
	ncm_user user_;

	//Subject we're current imaging
	ncm_subject subject_;
	
	//Study we're imaging for
	QString study_name_;

	//: Index of Finite State Machine that defines current activity
  activity_state_enum activity_state_;

	//Database object
	QSqlDatabase database_;

	//Flag confirming connection to database
	bool database_connected_;

	//Icons for starting/pausing recording (we need to switch between these)
	QIcon start_icon_;
	QIcon pause_icon_;

	//----------------------------------------------------
	//Reviewing

	//:QSql models for selecting sessions/sequences to reviews
	QSqlQueryModel *sessionModel_;
	QSortFilterProxyModel *sessionFilterModel_;
	QSqlQueryModel *sequenceModel_;
	QSortFilterProxyModel *sequenceFilterModel_;

	//: Class to handle image I/O and processing
  ncm_qimagehandler imageProcessor_;

	//: annotation class that stores image properties
  ncm_annotation markup_;

  //: canvas on which to draw stuff
  //  Note that this appears after markup_ in the declaration since scene_ needs
  //  the data in markup_ to destruct itself and must therefore be destroyed
  //  first (i.e. declared later than markup_)
  ncm_qscene review_scene_;
	QGraphicsSimpleTextItem* no_mosaic_msg_;
};

#endif // NCM_QCAPTURE_GUI_H
